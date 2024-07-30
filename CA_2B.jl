using Tensors, ForwardDiff, Waveforms
using MaterialModelsBase
using Newton
import CairoMakie as CM

include("taskb1.jl")
include("taskb2.jl")


main_yield_surface()   
main_demo_chabocheHill()  