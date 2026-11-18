using JunctionAnalysis
using GLMakie

# Choose :J to visualise supercurrent J_x(y) or :A to visualise vector potential A_x(y)
quantity = :J

# Narrow-junction regime matching analytic J_x(y) = J_DJ * cos(d * B * y):
# d = thick + 2*kappa = 7.5 < length = 16, width = 12 > 2*kappa = 6, lambda_J ≈ 9.4
junc = get_junction(0.015, 0.01, 2^-3, barrier = Barrier(-1.0, 1.0, 5.0, 1.5, 12, 16, 3.0))
fig, buffer = plot_junction(junc; quantity = quantity)

display(fig)

while junc.mode != JcDone() && isopen(fig.scene)
    step_animation!(junc, buffer)
end

save("stationary_state.png", fig)
