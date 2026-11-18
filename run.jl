using JunctionAnalysis

# Choose :J to visualise supercurrent J_x(y) or :A to visualise vector potential A_x(y)
quantity = :J

junc = get_junction(1e-4, 1e-4, 2^-5, barrier = Barrier(-1.0, 1.0, 5.0, 0.5, 2, 6, 30))
fig, buffer = plot_junction(junc; quantity = quantity)

display(fig)

while junc.mode != JcDone() && isopen(fig.scene)
    step_animation!(junc, buffer)
end