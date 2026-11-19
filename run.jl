using JunctionAnalysis
using GLMakie

# narrow = 1, wide = 32
width = 32

# Choose :J to visualise supercurrent J_x(y) or :A to visualise vector potential A_x(y)
quantity = :J

junc = get_junction(0.1, 0.001, 2^-3, barrier = Barrier(-1.0, 1.0, 1.0, 1.0, width, 16, 3.0))
println("λⱼ = ", josephson_penetration_depth(junc))
fig, buffer = plot_junction(junc; quantity = quantity)

display(fig)

while junc.mode != JcDone() && isopen(fig.scene)
    step_animation!(junc, buffer)
end

save("stationary_state.png", fig)
