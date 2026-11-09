struct PotentialBuffer
    A_obs::Observable{Vector{Vec2f}}
    J_obs::Observable{Vector{Vec2f}}
    ϕ_obs::Observable{Matrix{Float64}}
    ψ_obs::Observable{Matrix{Float64}}
    A_buff::Vector{Vec2f}
    J_buff::Vector{Vec2f}
    ϕ_buff::Matrix{Float64}
    ψ_buff::Matrix{Float64}
    vec_dims::NTuple{2, Int}
    step::Integer
end

"copy the vector potential data into a buffer, formatted as a vector of Vec2f"
function update_A_data!(sim::Solver, buffer::PotentialBuffer)
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
    update_A_data!(junction.solver, buffer)
    update_J_data!(junction, buffer)
    update_ϕ_data!(junction.solver, buffer)
    update_ψ_data!(junction.solver, buffer)

    notify(buffer.A_obs)
    notify(buffer.J_obs)
    notify(buffer.ϕ_obs)
    notify(buffer.ψ_obs)
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
function plot_junction(junction::Junction)
    sim = junction.solver

    fig = Figure(size = (800, 800))
    ax1 = Axis(fig[1,1])

    dims = min_dims(state(junction).a.extent)
    step = arrow_step(dims)
    vector_grid = grid_points(dims, step)

    A_grid = empty_vec_points(dims, step)
    A_obs = Observable(A_grid)

    ϕ_grid = adapt(CPU(), copy(retrieve_ϕ(sim)))
    ϕ_obs = Observable(ϕ_grid)

    ψ_grid = adapt(CPU(), copy(retrieve_ψ(sim)))
    ψ_obs = Observable(ψ_grid)

    J_grid = copy(A_grid)
    J_obs = Observable(J_grid)

    buffer = PotentialBuffer(A_obs, J_obs, ϕ_obs, ψ_obs, A_grid, J_grid, ϕ_grid, ψ_grid, dims, step)

    heatmap!(ax1, buffer.ψ_obs, colorrange = (0.0, 1.0))
    arrows2d!(ax1, vector_grid, rescale_arrows!(J_obs), align = :center)
    
    hidedecorations!(ax1)
    return fig, buffer 
end

function step_animation!(junction::Junction, buffer::PotentialBuffer)
    step!(junction)
    update_observables!(junction, buffer)
end