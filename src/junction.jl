const Solver = MulTDGL.ImplicitLondonMultigridSolver{2, Float64, 
               MulTDGL.RectPrimalForm0Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm1Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm0Data{2, ComplexF64, CUDACore.CuArray{ComplexF64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm1Data{2, ComplexF64, CUDACore.CuArray{ComplexF64, 1, CUDACore.DeviceMemory}},
               MulTDGL.RectPrimalForm2Data2{Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}}}

const Field = MulTDGL.RectPrimalForm1Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}}
const SimState = MulTDGL.State{
    MulTDGL.RectPrimalForm0Data{2, ComplexF64, CUDACore.CuArray{ComplexF64, 1, CUDACore.DeviceMemory}},
    MulTDGL.RectPrimalForm1Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}},
    MulTDGL.RectPrimalForm0Data{2, Float64, CUDACore.CuArray{Float64, 1, CUDACore.DeviceMemory}}
}

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
    min_e::Float64
    last_stable_state::SimState
    last_stable_current::Float64
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

function save_stable_state!(junction::Junction)
    copy!(junction.last_stable_state, state(junction))
    junction.last_stable_current = junction.current
end

function restore_stable_state!(junction::Junction)
    copy!(state(junction), junction.last_stable_state)
    junction.current = junction.last_stable_current
    MulTDGL.reconfigure!(system(junction), state(junction))
    apply_bcs!(junction)
    calculate_supercurrent!(junction)
end

function step!(junction::Junction)
    Ec = 5e-6
    apply_bcs!(junction)
    step_data = MulTDGL.step!(junction.solver, junction.boundary, (junction.current, 0.0))
    calculate_supercurrent!(junction)
    abs_e = abs(step_data.e[1])
    junction.holdsteps += 1

    if junction.mode == JcInitHold()
        if junction.holdsteps > 1_000 || (junction.holdsteps >= 100 && abs_e < Ec)
            save_stable_state!(junction)
            junction.holdsteps = 0
            junction.min_e = Inf
            junction.mode = JcJHold()
        end
    elseif junction.mode == JcJHold()
        if junction.holdsteps >= 5
            junction.min_e = min(junction.min_e, abs_e)
        end

        if junction.holdsteps >= 10 && abs_e < Ec
            save_stable_state!(junction)
            junction.holdsteps = 0
            junction.min_e = Inf
            junction.current *= 1.025
            println(junction.current)
        elseif (junction.holdsteps >= 15 && abs_e > 1.05 * junction.min_e) || junction.holdsteps >= 3_000
            restore_stable_state!(junction)
            calculate_supercurrent!(junction)

            junction.holdsteps = 0
            junction.mode = JcDone()
            println("Reached critical current J_c = $(junction.current) (restored stationary state)")
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

    last_stable_state = copy(MulTDGL.state(solver))

    return Junction(barrier, b_field, current, solver, boundary, j_super, JcInitHold(), 0, Inf, last_stable_state, Float64(current))
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

function barrier_x_coords(junction::Junction)
    msh = mesh(junction)
    hx = msh.h[1]
    nx = elemextent(msh, (), 1)
    x_mid = barrier_center_x(junction)
    return [(i - x_mid) * hx for i in 1:nx]
end

function effective_thickness(barrier::Barrier)
    return barrier.thick + 2.0 * barrier.kappa
end

function josephson_penetration_depth(barrier::Barrier, J_DJ::Real)
    d = effective_thickness(barrier)
    return barrier.kappa / sqrt(d * J_DJ)
end

josephson_penetration_depth(j::Junction) = josephson_penetration_depth(j.barrier, junction_depairing_current(j))

function analytic_narrow_J(junction::Junction, y::AbstractVector{<:Real}; J_DJ::Real = junction_depairing_current(junction))
    d = effective_thickness(junction.barrier)
    k = d * junction.field
    return J_DJ .* cos.(k .* y)
end

function analytic_narrow_A(junction::Junction, y::AbstractVector{<:Real})
    return -junction.field .* y
end

function analytic_wide_Jc(junction::Junction; J_DJ::Real = junction_depairing_current(junction))
    barrier = junction.barrier
    d = effective_thickness(barrier)
    w = barrier.width
    κ = barrier.kappa
    μ0 = 1.0 / κ^2
    λ_J = josephson_penetration_depth(barrier, J_DJ)
    B_0 = 2.0 / (d * λ_J)
    return max(0.0, (2.0 / (w * μ0)) * (B_0 - abs(junction.field)))
end

function wide_phase_profile(junction::Junction, _, y_point::Real, J_DJ::Real)
    barrier = junction.barrier
    d = effective_thickness(barrier)
    w = barrier.width
    κ = barrier.kappa
    μ0 = 1.0 / κ^2
    λ_J = josephson_penetration_depth(barrier, J_DJ)
    B_0 = 2.0 / (d * λ_J)

    B_app = junction.field
    J = analytic_wide_Jc(junction; J_DJ = J_DJ)
    db_half = 0.5 * μ0 * J * w

    B_plus = B_app + db_half
    B_minus = B_app - db_half

    # assumes Meissner state and that B can't exceed B_0
    ϕ_plus = 2.0 * asin(clamp(B_plus / B_0, -1.0, 1.0))
    ϕ_minus = 2.0 * asin(clamp(-B_minus / B_0, -1.0, 1.0))

    # tan and atan come from integrating dΔϕ/dy = 2/λⱼ * sin(Δϕ/2)
    t_plus = tan(ϕ_plus / 4.0)
    t_minus = tan(ϕ_minus / 4.0)

    Δϕ_plus = 4.0 * atan(t_plus * exp(-(w / 2.0 - y_point) / λ_J))
    Δϕ_minus = 4.0 * atan(t_minus * exp(-(y_point + w / 2.0) / λ_J))
    return Δϕ_plus + Δϕ_minus
end

function wide_field_at_y(junction::Junction, y_point::Real, J_DJ::Real = junction_depairing_current(junction))
    barrier = junction.barrier
    d = effective_thickness(barrier)
    w = barrier.width
    κ = barrier.kappa
    μ0 = 1.0 / κ^2
    λ_J = josephson_penetration_depth(barrier, J_DJ)
    B_0 = 2.0 / (d * λ_J)

    B_app = junction.field
    J = analytic_wide_Jc(junction; J_DJ = J_DJ)
    db_half = 0.5 * μ0 * J * w

    B_plus = B_app + db_half
    B_minus = B_app - db_half

    ϕ_plus = 2.0 * asin(clamp(B_plus / B_0, -1.0, 1.0))
    ϕ_minus = 2.0 * asin(clamp(-B_minus / B_0, -1.0, 1.0))

    t_plus = tan(ϕ_plus / 4.0)
    t_minus = tan(ϕ_minus / 4.0)

    Δϕ_plus = 4.0 * atan(t_plus * exp(-(w / 2.0 - y_point) / λ_J))
    Δϕ_minus = 4.0 * atan(t_minus * exp(-(y_point + w / 2.0) / λ_J))

    return B_0 * (sin(Δϕ_plus / 2.0) - sin(Δϕ_minus / 2.0))
end

function transverse_Bz_profile(junction::Junction, x::AbstractVector{<:Real}, B_center::Real)
    half_t = 0.5 * junction.barrier.thick
    κ = junction.barrier.kappa
    return [abs(xi) <= half_t ? Float64(B_center) : Float64(B_center) * exp(-(abs(xi) - half_t) / κ) for xi in x]
end

function analytic_narrow_Bz(junction::Junction, x::AbstractVector{<:Real})
    return transverse_Bz_profile(junction, x, junction.field)
end

function analytic_wide_Bz(junction::Junction, x::AbstractVector{<:Real}; J_DJ::Real = junction_depairing_current(junction))
    B_center = wide_field_at_y(junction, 0.0, J_DJ)
    return transverse_Bz_profile(junction, x, B_center)
end

function _electrode_f_infinity_sq(j_tilde_sq::Real)
    j_tilde_sq <= 0.0 && return 1.0
    max_j_sq = 4.0 / 27.0
    j_tilde_sq >= max_j_sq && return 2.0 / 3.0
    θ = acos(clamp(1.0 - 27.0 * j_tilde_sq / 2.0, -1.0, 1.0))
    return (1.0 + 2.0 * cos(θ / 3.0)) / 3.0
end

function _barrier_profile_and_deriv(u::Real, f0_sq::Real, j_u_sq::Real)
    f = sqrt(max(f0_sq, 1e-30))
    abs(u) < 1e-12 && return (f * f, 0.0)
    n_steps = max(32, ceil(Int, abs(u) / 0.005))
    du = abs(u) / n_steps
    df = 0.0

    @inline function accel(f_val)
        f2 = max(f_val * f_val, 1e-30)
        return (1.0 + f2 - j_u_sq / (f2 * f2)) * f_val
    end

    for _ in 1:n_steps
        k1_f = df
        k1_df = accel(f)

        f2 = f + 0.5 * du * k1_f
        df2 = df + 0.5 * du * k1_df
        k2_f = df2
        k2_df = accel(f2)

        f3 = f + 0.5 * du * k2_f
        df3 = df + 0.5 * du * k2_df
        k3_f = df3
        k3_df = accel(f3)

        f4 = f + du * k3_f
        df4 = df + du * k3_df
        k4_f = df4
        k4_df = accel(f4)

        f += (du / 6.0) * (k1_f + 2.0 * k2_f + 2.0 * k3_f + k4_f)
        df += (du / 6.0) * (k1_df + 2.0 * k2_df + 2.0 * k3_df + k4_df)
    end
    return (f * f, df)
end

function _barrier_profile(u::Real, f0_sq::Real, j_u_sq::Real)
    return _barrier_profile_and_deriv(u, f0_sq, j_u_sq)[1]
end

function _barrier_boundary_residual(f0_sq::Real, u_half::Real, j_u_sq::Real, v_sq::Real, s::Real, f_inf_sq::Real, j_tilde_sq::Real)
    f_d2_sq, df_d2 = _barrier_profile_and_deriv(u_half, f0_sq, j_u_sq)
    f_s_d2_sq = f_d2_sq / (-s)
    f_s_d2_sq > f_inf_sq && return Inf
    lhs = (df_d2 * df_d2) / (-s)
    rhs = 0.5 * v_sq * (2.0 * (f_inf_sq - f_s_d2_sq) - (f_inf_sq^2 - f_s_d2_sq^2) - 2.0 * j_tilde_sq * (1.0 / max(f_s_d2_sq, 1e-30) - 1.0 / f_inf_sq))
    return lhs - rhs
end

function _min_f0_sq(j_u_sq::Real, f_max_sq::Real)
    j_u_sq <= 1e-30 && return 0.0
    lo_a = 0.0
    lo_b = f_max_sq
    for _ in 1:60
        mid = 0.5 * (lo_a + lo_b)
        if 1.0 + mid - j_u_sq / max(mid * mid, 1e-30) < 0.0
            lo_a = mid
        else
            lo_b = mid
        end
    end
    return lo_b
end

function _solve_f0_sq(u_half::Real, j_u_sq::Real, v_sq::Real, s::Real, f_inf_sq::Real, j_tilde_sq::Real)
    f_max_sq = (-s) * f_inf_sq
    lo = _min_f0_sq(j_u_sq, f_max_sq)

    hi_a = lo
    hi_b = f_max_sq
    if _barrier_profile(u_half, hi_b, j_u_sq) / (-s) > f_inf_sq
        for _ in 1:60
            mid = 0.5 * (hi_a + hi_b)
            if _barrier_profile(u_half, mid, j_u_sq) / (-s) > f_inf_sq
                hi_b = mid
            else
                hi_a = mid
            end
        end
        hi = hi_a
    else
        hi = hi_b
    end

    lo >= hi && return lo

    r_lo = _barrier_boundary_residual(lo, u_half, j_u_sq, v_sq, s, f_inf_sq, j_tilde_sq)
    (isfinite(r_lo) && r_lo >= 0.0) && return lo

    bracket_lo = lo
    bracket_hi = hi
    for _ in 1:60
        mid = 0.5 * (bracket_lo + bracket_hi)
        r_mid = _barrier_boundary_residual(mid, u_half, j_u_sq, v_sq, s, f_inf_sq, j_tilde_sq)
        if !isfinite(r_mid) || r_mid > 0.0
            bracket_hi = mid
        else
            bracket_lo = mid
        end
    end
    return 0.5 * (bracket_lo + bracket_hi)
end

function _electrode_profile(x_tilde_dist::Real, f_s_d2_sq::Real, f_inf_sq::Real, j_tilde_sq::Real)
    x_tilde_dist <= 1e-12 && return f_s_d2_sq
    f_s_d2_sq >= f_inf_sq && return f_inf_sq

    n_steps = max(64, ceil(Int, x_tilde_dist / 0.01))
    dx = x_tilde_dist / n_steps
    f = sqrt(max(f_s_d2_sq, 1e-30))
    rhs_0 = 0.5 * (2.0 * (f_inf_sq - f_s_d2_sq) - (f_inf_sq^2 - f_s_d2_sq^2) - 2.0 * j_tilde_sq * (1.0 / max(f_s_d2_sq, 1e-30) - 1.0 / f_inf_sq))
    df = sqrt(max(0.0, rhs_0))

    @inline function accel(f_val)
        f2 = clamp(f_val * f_val, 1e-30, f_inf_sq)
        return -(1.0 - f2 - j_tilde_sq / (f2 * f2)) * f_val
    end

    for _ in 1:n_steps
        k1_f = df
        k1_df = accel(f)

        f2 = f + 0.5 * dx * k1_f
        df2 = max(0.0, df + 0.5 * dx * k1_df)
        k2_f = df2
        k2_df = accel(f2)

        f3 = f + 0.5 * dx * k2_f
        df3 = max(0.0, df + 0.5 * dx * k2_df)
        k3_f = df3
        k3_df = accel(f3)

        f4 = f + dx * k3_f
        df4 = max(0.0, df + dx * k3_df)
        k4_f = df4
        k4_df = accel(f4)

        f = min(sqrt(f_inf_sq), f + (dx / 6.0) * (k1_f + 2.0 * k2_f + 2.0 * k3_f + k4_f))
        df = max(0.0, df + (dx / 6.0) * (k1_df + 2.0 * k2_df + 2.0 * k3_df + k4_df))
    end
    return min(f_inf_sq, f * f)
end

function _has_valid_barrier_solution(j_eff::Real, one_minus_B::Real, bar::Barrier, xi_n::Real, u_half::Real, v_sq::Real, s::Real)
    j_tilde = j_eff / (one_minus_B^1.5)
    j_tilde_sq = j_tilde * j_tilde
    j_u = j_eff * bar.beta * sqrt(bar.mass) * (xi_n^3)
    j_u_sq = j_u * j_u

    f_inf_sq = _electrode_f_infinity_sq(j_tilde_sq)
    f_max_sq = (-s) * f_inf_sq
    lo = _min_f0_sq(j_u_sq, f_max_sq)
    r_lo = _barrier_boundary_residual(lo, u_half, j_u_sq, v_sq, s, f_inf_sq, j_tilde_sq)
    return isfinite(r_lo) && r_lo <= 0.0
end

function transverse_ψ_profile(junction::Junction, x::AbstractVector{<:Real}, J_local::Real, B_local::Real)
    bar = junction.barrier
    B = clamp(abs(B_local), 0.0, 1.0 - 1e-9)
    one_minus_B = 1.0 - B

    inv_xi_n_sq = -bar.mass * bar.alpha + B
    xi_n = sqrt(1.0 / inv_xi_n_sq)
    s = bar.beta * one_minus_B / (bar.alpha - B / bar.mass)
    v = bar.mass * xi_n * sqrt(one_minus_B)
    v_sq = v * v

    half_t = 0.5 * bar.thick
    u_half = half_t / xi_n
    sqrt_one_minus_B = sqrt(one_minus_B)

    J_DJ = junction_depairing_current(junction)
    J_eff = clamp(abs(J_local), 0.0, J_DJ)

    if !_has_valid_barrier_solution(J_eff, one_minus_B, bar, xi_n, u_half, v_sq, s)
        j_lo = 0.0
        j_hi = J_eff
        for _ in 1:50
            j_mid = 0.5 * (j_lo + j_hi)
            if _has_valid_barrier_solution(j_mid, one_minus_B, bar, xi_n, u_half, v_sq, s)
                j_lo = j_mid
            else
                j_hi = j_mid
            end
        end
        J_eff = j_lo
    end

    j_tilde = J_eff / (one_minus_B^1.5)
    j_tilde_sq = j_tilde * j_tilde

    j_u = J_eff * bar.beta * sqrt(bar.mass) * (xi_n^3)
    j_u_sq = j_u * j_u

    f_inf_sq = _electrode_f_infinity_sq(j_tilde_sq)
    f0_sq = _solve_f0_sq(u_half, j_u_sq, v_sq, s, f_inf_sq, j_tilde_sq)
    f_d2_sq = _barrier_profile(u_half, f0_sq, j_u_sq)
    f_s_d2_sq = min(f_inf_sq, f_d2_sq / (-s))

    scale_n = (-bar.alpha + B / bar.mass) / bar.beta
    scale_s = one_minus_B

    out = Vector{Float64}(undef, length(x))
    for (idx, xi) in pairs(x)
        ax = abs(xi)
        if ax <= half_t
            u = ax / xi_n
            f_n_sq = _barrier_profile(u, f0_sq, j_u_sq)
            out[idx] = scale_n * f_n_sq
        else
            x_tilde_dist = (ax - half_t) * sqrt_one_minus_B
            f_s_sq = _electrode_profile(x_tilde_dist, f_s_d2_sq, f_inf_sq, j_tilde_sq)
            out[idx] = scale_s * f_s_sq
        end
    end
    return out
end

function analytic_narrow_ψ(junction::Junction, x::AbstractVector{<:Real}; J_DJ::Real = junction_depairing_current(junction))
    return transverse_ψ_profile(junction, x, J_DJ, junction.field)
end

function analytic_wide_ψ(junction::Junction, x::AbstractVector{<:Real}; J_DJ::Real = junction_depairing_current(junction))
    Δϕ_0 = wide_phase_profile(junction, nothing, 0.0, J_DJ)
    J_0 = J_DJ * sin(Δϕ_0)
    B_0 = wide_field_at_y(junction, 0.0, J_DJ)
    return transverse_ψ_profile(junction, x, J_0, B_0)
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
