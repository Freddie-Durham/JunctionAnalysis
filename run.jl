using JunctionAnalysis

junc = get_junction(0.001, 1e-4, barrier = Barrier(-1.0, 1.0, 1.0, 4, 256, 256, 30))
fig, buffer = plot_junction(junc)

display(fig)

while junc.mode != JcDone() && isopen(fig.scene)
    step_animation!(junc, buffer)
end