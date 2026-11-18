struct PotentialBuffer
    A_obs::Observable{Vector{Vec2f}}
    J_obs::Observable{Vector{Vec2f}}
    Bz_obs::Observable{Matrix{Float64}}
    ϕ_obs::Observable{Matrix{Float64}}
    ψ_obs::Observable{Matrix{Float64}}
    sim_1d_obs::Observable{Vector{Float64}}
    narrow_1d_obs::Observable{Vector{Float64}}
    wide_1d_obs::Observable{Vector{Float64}}
    A_buff::Vector{Vec2f}
    J_buff::Vector{Vec2f}
    Bz_buff::Matrix{Float64}
    ϕ_buff::Matrix{Float64}
    ψ_buff::Matrix{Float64}
    sim_1d_buff::Vector{Float64}
    narrow_1d_buff::Vector{Float64}
    wide_1d_buff::Vector{Float64}
    y_coords::Vector{Float64}
    x_mid::Int
    quantity::Symbol
    vec_dims::NTuple{2, Int}
    step::Integer
end

"copy the vector potential data into a buffer, formatted as a vector of Vec2f"
function update_A_data!(junction::Junction, buffer::PotentialBuffer)
    sim = junction.solver
    A_edges = adapt(CPU(), MulTDGL.state(sim).a)
    step = buffer.step
    dims = buffer.vec_dims
   
    idx = 1
    for x in 1:step:dims[1], y in 1:step:dims[2]
        i = CartesianIndex(x, y)
        ex = MulTDGL.Primal(MulTDGL.RectIndex(i, (1,)))
        ey = MulTDGL.Primal(MulTDGL.RectIndex(i, (2,)))

        buffer.A_buff[idx] = Vec2f(A_edges[ex], A_edges[ey])
        idx += 1
    end

    hx, hy = mesh(junction).h
    inv_area = 1.0 / (hx * hy)
    A_x = MulTDGL.asarray(A_edges, 1)
    A_y = MulTDGL.asarray(A_edges, 2)
    nx, ny = size(buffer.Bz_buff)
    @inbounds for y in 1:ny, x in 1:nx
        x_next = x == nx ? 1 : x + 1
        buffer.Bz_buff[x, y] = ((A_y[x_next, y] - A_y[x, y]) - (A_x[x, y + 1] - A_x[x, y])) * inv_area
    end

    if buffer.quantity == :A
        buffer.sim_1d_buff .= @view(A_x[buffer.x_mid, :]) ./ hx
        buffer.wide_1d_buff .= analytic_wide_A(junction, buffer.y_coords)
    end
end

function update_J_data!(junction::Junction, buffer::PotentialBuffer)
    J_edges = adapt(CPU(), junction.supercurrent)
    step = buffer.step
    dims = buffer.vec_dims
   
    idx = 1
    for x in 1:step:dims[1], y in 1:step:dims[2]
        i = CartesianIndex(x, y)
        ex = MulTDGL.Primal(MulTDGL.RectIndex(i, (1,)))
        ey = MulTDGL.Primal(MulTDGL.RectIndex(i, (2,)))

        buffer.J_buff[idx] = Vec2f(J_edges[ex], J_edges[ey])
        idx += 1
    end

    if buffer.quantity == :J
        hx = mesh(junction).h[1]
        J_x = MulTDGL.asarray(J_edges, 1)
        buffer.sim_1d_buff .= @view(J_x[buffer.x_mid, :]) ./ hx
        buffer.wide_1d_buff .= analytic_wide_J(junction, buffer.y_coords)
    end
end

retrieve_ϕ(sim) = asarray(MulTDGL.state(sim).φ)

"copy the scalar potential data into a buffer"
function update_ϕ_data!(sim::Solver, buffer::PotentialBuffer)
    ϕ_data = retrieve_ϕ(sim)
    copyto!(buffer.ϕ_buff, ϕ_data)
end

retrieve_ψ(sim) = abs.(asarray(MulTDGL.state(sim).ψ)) .^2

"copy the order parameter data into a buffer"
function update_ψ_data!(sim::Solver, buffer::PotentialBuffer)
    ψ_data = retrieve_ψ(sim)
    copyto!(buffer.ψ_buff, ψ_data)
end

"update buffers with latest data and use notify to update Observables"
function update_observables!(junction::Junction, buffer::PotentialBuffer)
    update_A_data!(junction, buffer)
    update_J_data!(junction, buffer)
    update_ϕ_data!(junction.solver, buffer)
    update_ψ_data!(junction.solver, buffer)

    notify(buffer.A_obs)
    notify(buffer.J_obs)
    notify(buffer.Bz_obs)
    notify(buffer.ϕ_obs)
    notify(buffer.ψ_obs)
    notify(buffer.sim_1d_obs)
    notify(buffer.wide_1d_obs)
end

const ARROW_COUNT = 16
arrow_step(dims) = ceil(Int, dims[1] / ARROW_COUNT)
min_dims(dims) = reduce((cur, next) -> min.(cur, next), dims)

grid_points(dims, step) = vec([Point2f(x + step / 2, y + step / 2) for x in 1:step:dims[1], y in 1:step:dims[2]])
empty_vec_points(dims, step) = vec([Vec2f(0, 0) for _ in 1:step:dims[1], _ in 1:step:dims[2]])

function rescale_arrows!(edge_obs)
    return @lift(begin
        vectors = $edge_obs
        max_mag = maximum(v -> sqrt(v[1]^2 + v[2]^2), vectors)
        
        if max_mag > 0.0
            (vectors ./ max_mag) .* 10.0
        else
            vectors
        end
    end)
end

"initialise heatmap of 2D node data and arrows! plot of edge data that can be updated using Observables"
function plot_junction(junction::Junction; quantity::Symbol = :J)
    quantity in (:J, :A) || throw(ArgumentError("quantity must be :J or :A, got $(repr(quantity))"))
    sim = junction.solver

    fig = Figure(size = (1400, 700))
    ax1 = Axis(fig[1, 1], aspect = DataAspect())

    dims = min_dims(state(junction).a.extent)
    step = arrow_step(dims)

    A_grid = empty_vec_points(dims, step)
    A_obs = Observable(A_grid)

    face_dims = elemextent(mesh(junction), (1, 2))
    Bz_grid = zeros(Float64, face_dims)
    Bz_obs = Observable(Bz_grid)

    ϕ_grid = adapt(CPU(), copy(retrieve_ϕ(sim)))
    ϕ_obs = Observable(ϕ_grid)

    ψ_grid = adapt(CPU(), copy(retrieve_ψ(sim)))
    ψ_obs = Observable(ψ_grid)

    J_grid = copy(A_grid)
    J_obs = Observable(J_grid)

    y_coords = barrier_y_coords(junction)
    x_mid = barrier_center_x(junction)

    sim_1d_buff = zeros(Float64, length(y_coords))
    J_DJ = junction_depairing_current(junction)
    narrow_1d_buff = quantity == :J ? analytic_narrow_J(junction, y_coords; J_DJ = J_DJ) : analytic_narrow_A(junction, y_coords)
    wide_1d_buff = quantity == :J ? analytic_wide_J(junction, y_coords; J_DJ = J_DJ) : analytic_wide_A(junction, y_coords; J_DJ = J_DJ)

    sim_1d_obs = Observable(sim_1d_buff)
    narrow_1d_obs = Observable(narrow_1d_buff)
    wide_1d_obs = Observable(wide_1d_buff)

    buffer = PotentialBuffer(
        A_obs, J_obs, Bz_obs, ϕ_obs, ψ_obs,
        sim_1d_obs, narrow_1d_obs, wide_1d_obs,
        A_grid, J_grid, Bz_grid, ϕ_grid, ψ_grid,
        sim_1d_buff, narrow_1d_buff, wide_1d_buff,
        y_coords, x_mid, quantity, dims, step
    )

    heatmap!(ax1, buffer.ψ_obs, colorrange = (0.0, 1.0))
    nx, ny_verts = size(ψ_grid)
    contour!(
        ax1,
        range(1, nx, length = face_dims[1]),
        range(1, ny_verts, length = face_dims[2]),
        buffer.Bz_obs;
        levels = 15,
        color = :white,
        linewidth = 1.5
    )
    hidedecorations!(ax1)

    ylabel_str = quantity == :J ? "J_x" : "A_x"
    title_str = quantity == :J ? "Barrier Supercurrent Profile J_x(y)" : "Barrier Vector Potential Profile A_x(y)"
    ax2 = Axis(fig[1, 2], xlabel = "y (centered on barrier)", ylabel = ylabel_str, title = title_str)

    lines!(ax2, y_coords, sim_1d_obs, label = "Simulated", linewidth = 2)
    lines!(ax2, y_coords, narrow_1d_obs, label = "Analytic (Narrow)", linestyle = :dash, linewidth = 2)
    lines!(ax2, y_coords, wide_1d_obs, label = "Analytic (Wide)", linestyle = :dot, linewidth = 2)
    axislegend(ax2, position = :rt)

    return fig, buffer 
end

function step_animation!(junction::Junction, buffer::PotentialBuffer)
    step!(junction)
    update_observables!(junction, buffer)
end