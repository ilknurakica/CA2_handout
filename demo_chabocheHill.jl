using Tensors, ForwardDiff, Waveforms
using MaterialModelsBase
using Newton
import CairoMakie as CM

# Include the supplied files
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "perfect_plasticity.jl"))
include(joinpath(@__DIR__, "chaboche_hill.jl"))
include(joinpath(@__DIR__, "uniaxial_stress.jl"))

# Helper function to get a triangle wave strain time history
function get_cycles(;amplitude=1.0, num_cycles=1, num_steps_per_cycle=100)
    return amplitude*trianglewave1.(range(0.0, num_cycles, num_steps_per_cycle*num_cycles))
end

# Define material

G_Hill = 0.125
H_Hill = 0.625
F_Hill = 0.375
K_Hill = 1.500
L_Hill = 1.500
M_Hill = 1.725
"""
F_Hill=0.5
G_Hill=0.5
H_Hill=0.5
K_Hill=3.0
L_Hill=3.0
M_Hill=3.0
"""
m_chaboche_hill = ChabocheHill(G=80e3, K=160e3, Y=100.0, Hiso=0.0, κ∞=400.0, Hkin=0.0, β∞=50.0,
                               G_Hill=G_Hill, H_Hill=H_Hill, F_Hill=F_Hill,
                               K_Hill=K_Hill, L_Hill=L_Hill, M_Hill=M_Hill)
#Hiso=100e3, κ∞=400.0, Hkin=400e3, β∞=50.0
# Function to plot results
function plot_results(m::AbstractMaterial, ϵ11::Vector, ϵ12::Vector; kwargs...)
    fig = CM.Figure()
    ax1 = CM.Axis(fig[1,1]; xlabel = "Time", ylabel = "σ₁₁ [MPa]", title="σ₁₁ vs Time")
    ax2 = CM.Axis(fig[1,2]; xlabel = "Time", ylabel = "σ₁₂ [MPa]", title="σ₁₂ vs Time")
    
plot_results!(ax1, ax2, m, ϵ11, ϵ12; label = "case 1")
    
    return fig, ax1, ax2
end

function plot_results!(ax1, ax2,  m::AbstractMaterial, ϵ11, ϵ12; label="", t = collect(range(0,1,length(ϵ11))))
    σ11, σ12 = uniaxial_stress(m, ϵ11, ϵ12, t)
  
    CM.lines!(ax1, t, σ11; label=label)
    CM.lines!(ax2, t, σ12; label=label)
    
    CM.axislegend(ax1; position=:rt)
    CM.axislegend(ax2; position=:rt)
    return nothing
end


# Case 1: ϵ11 ≠ 0, ϵ12 = 0
ϵ11_case1 = get_cycles(amplitude=1, num_cycles=1, num_steps_per_cycle=100)
ϵ12_case1 = zeros(100)

fig, ax1, ax2 = plot_results(m_chaboche_hill, ϵ11_case1, ϵ12_case1)



# Case 2: ϵ11 = ϵ12 ≠ 0
ϵ11_case2 = get_cycles(amplitude=1, num_cycles=1, num_steps_per_cycle=100)
ϵ12_case2 = ϵ11_case2

plot_results!(ax1, ax2, m_chaboche_hill, ϵ11_case2, ϵ12_case2; label="case 2")

# Case 3: ϵ11 = 0, ϵ12 ≠ 0
ϵ11_case3 = zeros(100)
ϵ12_case3 = get_cycles(amplitude=1, num_cycles=1, num_steps_per_cycle=100)

plot_results!(ax1, ax2, m_chaboche_hill, ϵ11_case3, ϵ12_case3; label="case 3")

CM.axislegend(ax; position=:rb)
fig


