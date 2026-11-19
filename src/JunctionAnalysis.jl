module JunctionAnalysis

using CUDA
using KernelAbstractions
using Adapt
using MulTDGL
using GLMakie
using Random

include("junction.jl")
include("visualise.jl")

export Barrier, get_junction, plot_junction, 
       state, mesh, step!, step_animation!,
       JcInitHold, JcJHold, JcDone,
       barrier_center_x, barrier_y_coords,
       junction_depairing_current, josephson_penetration_depth,
       analytic_narrow_J, analytic_wide_J,
       analytic_narrow_A, analytic_wide_A

end # module