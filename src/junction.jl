const Solver = MulTDGL.ImplicitLondonMultigridSolver{2, Float64, 
               MulTDGL.RectPrimalForm0Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm1Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm0Data{2, ComplexF64, CUDACore.CuArray{ComplexF64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm1Data{2, ComplexF64, CUDACore.CuArray{ComplexF64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm2Data2{Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}}}

const Field = MulTDGL.RectPrimalForm1Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}}

struct JcInitHold end
struct JcJHold end
struct JcDone end
const JcMode = Union{JcInitHold, JcJHold, JcDone}

linear_interp(min_val, max_val, val) = min_val + (max_val - min_val) * val

function interp_edgesY(mesh, periodic_y)
    dims = size(mesh)
    new_mesh = zeros((dims[1], dims[2] - 1 + periodic_y))
    for x in 1:dims[1]
        for y in 2:dims[2]
            new_mesh[x, y-1] = (mesh[x, y-1] + mesh[x, y]) / 2
        end
        if periodic_y
            new_mesh[x, end] = (mesh[x, end] + mesh[x, 1]) / 2
        end
    end
    return new_mesh
end

function interp_edgesX(mesh, periodic_x)
    dims = size(mesh)
    new_mesh = zeros((dims[1] - 1 + periodic_x, dims[2]))
    for y in 1:dims[2]
        for x in 2:dims[1]
            new_mesh[x-1, y] = (mesh[x-1, y] + mesh[x, y]) / 2
        end
        if periodic_x
            new_mesh[end, y] = (mesh[end, y] + mesh[1, y]) / 2
        end
    end
    return new_mesh
end

struct Barrier
    alpha::Float64
    beta::Float64
    mass::Float64
    thick::Float64
    width::Int32
    length::Int32
    kappa::Float64
end

Barrier() = Barrier(-1.0, 1.0, 1.0, 2, 40, 40, 30)

mutable struct Junction
    barrier::Barrier
    field::Float64
    current::Float64
    solver::Solver
    boundary::Field
    supercurrent::Field
    mode::JcMode
    holdsteps::Int64
end

state(junction::Junction) = MulTDGL.state(junction.solver)
mesh(junction::Junction) = MulTDGL.mesh(junction.solver)
parameters(junction::Junction) = MulTDGL.parameters(junction.solver)
system(junction::Junction) = MulTDGL.system(junction.solver)

function calculate_supercurrent!(junction::Junction)
    st = state(junction)
    msh = mesh(junction)
    sys = system(junction)

    jₛ!(junction.supercurrent, msh, sys.mat.m⁻¹, sys.u, st.ψ)
end

function apply_bcs!(junction::Junction)
    j = junction.current
    params = parameters(junction)
    msh = mesh(junction)

    if !msh.periodic[2]
        b_field = junction.field
        
        db = j * msh.extent[2] * msh.h[2] / params.κ^2
        MulTDGL.asarray(junction.boundary, 1)[:, 1] .= (b_field - db/2) * 2 * msh.h[1] / msh.h[2]
        MulTDGL.asarray(junction.boundary, 1)[:, end] .= -(b_field + db/2) * 2 * msh.h[1] / msh.h[2]

    else all(msh.periodic)
        MulTDGL.asarray(junction.boundary, 1) .= -j * msh.h[1] / params.κ^2
    end
end

function step!(junction::Junction)
    apply_bcs!(junction)
    step_data = MulTDGL.step!(junction.solver, junction.boundary, (junction.current, 0.0))
    calculate_supercurrent!(junction)
    e_field = step_data.e[1]
    junction.holdsteps += 1

    if junction.mode == JcInitHold()
        if junction.holdsteps > 1_000
            junction.holdsteps = 0
            junction.mode = JcJHold()
        end
    elseif junction.mode == JcJHold()
        if e_field < 1e-5
            junction.holdsteps = 0
            junction.current *= 1.025
            println(junction.current)
        elseif junction.holdsteps >= 10_000
            junction.holdsteps = 0
            junction.mode = JcDone()
        end
    end
end

function apply_pattern(h, barrier, veclen)
    init_grid = zeros(Float64, veclen)

    junc_thick = Int(barrier.thick / h)

    half_length = floor(Int, (barrier.length / (2 * h) - barrier.thick / h))

    init_grid[half_length+1:half_length+junc_thick, :] .= 1.0
    return init_grid
end

function get_mesh(h, barrier, periodic)
    extent = Int.((barrier.length, barrier.width) ./ h)
    periodic = (true, periodic)
    hs = (h, h)
    return MulTDGL.RectMesh(extent, periodic, hs)
end

function get_material(h, barrier, mesh, backend)
    weights = apply_pattern(h, barrier, elemextent(mesh, ()))

    start_α = map(w -> linear_interp(1.0, barrier.alpha, w), weights)
    start_β = map(w -> linear_interp(1.0, barrier.beta, w), weights)

    start_σX = KernelAbstractions.ones(backend, Float64, elemextent(mesh, (1,)))
    start_σY = KernelAbstractions.ones(backend, Float64, elemextent(mesh, (2,)))

    inv_mass = 1 / barrier.mass
    start_m⁻¹X = map(w -> linear_interp(1.0, inv_mass, w), interp_edgesX(weights, mesh.periodic[1]))
    start_m⁻¹Y = map(w -> linear_interp(1.0, inv_mass, w), interp_edgesY(weights, mesh.periodic[2]))

    α = RectPrimalForm0Data(mesh, adapt(backend, reshape(start_α, :)))
    β = RectPrimalForm0Data(mesh, adapt(backend, reshape(start_β, :)))
    σ = RectPrimalForm1Data(mesh, adapt(backend, vcat(reshape(start_σX, :), reshape(start_σY, :)))) 
    m⁻¹ = RectPrimalForm1Data(mesh, adapt(backend, vcat(reshape(start_m⁻¹X, :), reshape(start_m⁻¹Y, :)))) 

    return MulTDGL.Material(α, β, m⁻¹, σ)
end

function get_system(mesh, params, material, backend, b_field)
    n1 = elemcount(mesh, Val(1))
    n2 = elemcount(mesh, Val(2))

    c = RectPrimalForm1Data(mesh, KernelAbstractions.zeros(backend, Float64, n1))
    dc = RectPrimalForm2Data2(mesh, KernelAbstractions.zeros(backend, Float64, n2))
    u = RectPrimalForm1Data(mesh, KernelAbstractions.ones(backend, Complex{Float64}, n1))

    if mesh.periodic[2]
        set_fluxons!(c, dc, b_field, mesh)
    end

    return System(mesh, params, material, c, dc, u)
end

function get_state(mesh, backend)
    n0 = elemcount(mesh, Val(0))
    n1 = elemcount(mesh, Val(1))

    ψ = RectPrimalForm0Data(mesh, KernelAbstractions.ones(backend, Complex{Float64}, n0))
    a = RectPrimalForm1Data(mesh, KernelAbstractions.zeros(backend, Float64, n1))
    φ = RectPrimalForm0Data(mesh, KernelAbstractions.zeros(backend, Float64, n0))

    return MulTDGL.State(ψ, a, φ)
end

function get_boundary(solver)
    boundary = MulTDGL.similar(MulTDGL.state(solver).a)
    MulTDGL.data(boundary) .= zero(eltype(MulTDGL.data(boundary)))
    return boundary
end

function get_junction(b_field, current, h=0.5; barrier = Barrier()) 
    backend = CUDABackend()

    Δt = 0.5
    tol = 1e-3
    levelcount = 3

    params = MulTDGL.Parameters(Δt, barrier.kappa)
    mesh = get_mesh(h, barrier, false)
    material = get_material(h, barrier, mesh, backend)
    system = get_system(mesh, params, material, backend, b_field)
    state = get_state(mesh, backend)

    solver = MulTDGL.ImplicitLondonMultigridSolver(system, state, tol, levelcount)
    boundary = get_boundary(solver)

    j_super = MulTDGL.similar(MulTDGL.state(solver).a)
    MulTDGL.data(j_super) .= zero(eltype(MulTDGL.data(j_super)))

    return Junction(barrier, b_field, current, solver, boundary, j_super, JcInitHold(), 0)
end

const DEPAIRING_CURRENT_1D = 2 / (3 * sqrt(3))

function junction_depairing_current(junction::Junction)
    bar = junction.barrier
    B = junction.field

    xi = sqrt.(1.0 ./ (-bar.mass .* bar.alpha .+ B))
    s = bar.beta .* (1.0 .- B) ./ (bar.alpha .- B ./ bar.mass)
    v = bar.mass .* xi .* sqrt.(1.0 .- B)

    f_d2 = if isapprox(bar.mass, bar.beta)
        (1.0 .- B) ./ (2.0 .- 2.0 .* bar.alpha)
    else
        (v.^2 .+ 1.0 .- sqrt.(max.(0.0, 1.0 .+ (2.0 .- s) .* v.^2))) ./ (v.^2 .+ s)
    end

    fink = (1.0 .- sqrt.(max.(0.0, 1.0 .- s .* f_d2))) ./ (s .* v)
    boyd = (1.0 .- B).^(1.5)
    decay = exp.(-bar.thick ./ xi)

    return 4.0 .* fink .* boyd .* decay
end

function barrier_center_x(junction::Junction)
    h = mesh(junction).h[1]
    barrier = junction.barrier
    junc_thick = Int(barrier.thick / h)
    half_length = floor(Int, (barrier.length - 2 * barrier.thick) / (2 * h))
    return half_length + cld(junc_thick, 2)
end

function barrier_y_coords(junction::Junction)
    msh = mesh(junction)
    w = msh.extent[2] * msh.h[2]
    ny = elemextent(msh, (1,), 2)
    return collect(range(-w / 2, w / 2, length = ny))
end

function effective_thickness(barrier::Barrier)
    return Float64(barrier.thick) + 2.0 * barrier.kappa
end

function josephson_penetration_depth(barrier::Barrier, J_DJ::Real)
    d = effective_thickness(barrier)
    return barrier.kappa / sqrt(d * J_DJ)
end

function analytic_narrow_J(junction::Junction, y::AbstractVector{<:Real}; J_DJ::Real = junction_depairing_current(junction))
    d = effective_thickness(junction.barrier)
    k = d * junction.field
    return J_DJ .* cos.(k .* y)
end

function analytic_narrow_A(junction::Junction, y::AbstractVector{<:Real})
    return -junction.field .* y
end

function wide_phase_profile(junction::Junction, _, y_point::Real, J_DJ::Real)
    barrier = junction.barrier
    d = effective_thickness(barrier)
    w = Float64(barrier.width)
    κ = barrier.kappa
    μ0 = 1.0 / κ^2
    λ_J = josephson_penetration_depth(barrier, J_DJ)
    B_0 = 2.0 / (d * λ_J)

    B_app = junction.field
    J = junction.current
    db_half = 0.5 * μ0 * J * w

    B_plus = B_app + db_half
    B_minus = B_app - db_half

    ϕ_plus = 2.0 * asin(clamp(B_plus / B_0, -1.0, 1.0))
    ϕ_minus = 2.0 * asin(clamp(-B_minus / B_0, -1.0, 1.0))

    t_plus = tan(ϕ_plus / 4.0)
    t_minus = tan(ϕ_minus / 4.0)

    Δϕ_plus = 4.0 * atan(t_plus * exp(-(w / 2.0 - y_point) / λ_J))
    Δϕ_minus = 4.0 * atan(t_minus * exp(-(y_point + w / 2.0) / λ_J))
    return Δϕ_plus + Δϕ_minus
end

function wide_phase_profile(junction::Junction, y::AbstractVector{<:Real}, J_DJ::Real)
    return [wide_phase_profile(junction, y, yi, J_DJ) for yi in y]
end

function analytic_wide_J(junction::Junction, y::AbstractVector{<:Real}; J_DJ::Real = junction_depairing_current(junction))
    Δϕ = wide_phase_profile(junction, y, J_DJ)
    return J_DJ .* sin.(Δϕ)
end

function analytic_wide_A(junction::Junction, y::AbstractVector{<:Real}; J_DJ::Real = junction_depairing_current(junction))
    d = effective_thickness(junction.barrier)
    Δϕ = wide_phase_profile(junction, y, J_DJ)
    Δϕ_0 = wide_phase_profile(junction, y, 0.0, J_DJ)
    return -(Δϕ .- Δϕ_0) ./ d
end
