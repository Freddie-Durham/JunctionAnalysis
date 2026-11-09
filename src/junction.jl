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
    thick::Int32
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
    half_length = floor(Int, (barrier.length - 2 * barrier.thick) / (2 * h))
    coat_width = Int(5.0 / h)

    init_grid[half_length+1:half_length+junc_thick, :] .= 1.0
    init_grid[:, 1:coat_width] .= 1.0
    init_grid[:, end-coat_width+1:end] .= 1.0
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

    ψdata = Adapt.adapt(backend, exp.(Complex{Float64}(2π * im) .* rand(Xoshiro(1), Float64, n0)))
    ψ = RectPrimalForm0Data(mesh, ψdata)
    a = RectPrimalForm1Data(mesh, KernelAbstractions.zeros(backend, Float64, n1))
    φ = RectPrimalForm0Data(mesh, KernelAbstractions.zeros(backend, Float64, n0))

    return MulTDGL.State(ψ, a, φ)
end

function get_boundary(solver)
    boundary = MulTDGL.similar(MulTDGL.state(solver).a)
    MulTDGL.data(boundary) .= zero(eltype(MulTDGL.data(boundary)))
    return boundary
end

function get_junction(b_field, current; barrier = Barrier()) 
    backend = CUDABackend()

    h = 0.5
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

    return Junction(barrier, b_field, current, solver, boundary, j_super, JcInitHold(), 0)
end
