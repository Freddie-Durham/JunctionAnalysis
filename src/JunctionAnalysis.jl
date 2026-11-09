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
       JcInitHold, JcJHold, JcDone

end # module