using Tensors, ForwardDiff, Waveforms
using MaterialModelsBase
using Newton
import CairoMakie as CM
function taska4()
# Include the supplied files
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "perfect_plasticity.jl"))
include(joinpath(@__DIR__, "chaboche.jl"))
include(joinpath(@__DIR__, "uniaxial_stress.jl"))

# Helper function to get a triangle wave strain time history
function get_cycles(;amplitude=1.0, num_cycles=3, num_steps_per_cycle=100)
    return amplitude*trianglewave1.(range(0.0, num_cycles, num_steps_per_cycle*num_cycles))
end

# Define the material
m_perfect = PerfectPlasticity(G=80.e3, K=160.e3, Y=250.0)
m_chaboche = Chaboche(G=80.e3, K=160.e3, Y=250.0, Hiso=100e3, κ∞=400.0, Hkin=40e3, β∞=50.0)
# Hiso=0.0, κ∞=400.0, Hkin=0.0, β∞=50.0)

#Hiso=100e3, κ∞=400.0, Hkin=40e3, β∞=50.0)

# Define function to plot results
function plot_results(m, ϵ11; kwargs...)
    fig = CM.Figure();
    ax = CM.Axis(fig[1,1]; xlabel = "ϵ₁₁ [%]", ylabel = "σ₁₁ [MPa]")
plot_results!(ax, m, ϵ11; kwargs...)    
    return fig, ax
end
function plot_results!(ax, m::AbstractMaterial, ϵ11; label="", t = collect(range(0,1,length(ϵ11))))
    σ11 = uniaxial_stress_original(m, ϵ11, t)
    CM.lines!(ax, 100*ϵ11, σ11; label=label)
    return nothing
end

# Example how to plot 1 cycle, for different strain amplitudes and monotonic loading
function create_example_plot(m, percent)
    fig = CM.Figure((size = (1000, 600)));
    ax = CM.Axis(fig[1,1]; xlabel = "ϵ₁₁ [%]", ylabel = "σ₁₁ [MPa]")
    # fig, ax = plot_results(m, collect(range(0, 1.0/100, 200)); label="monotonic")
    for amplitude_percent in percent
        ϵ11 = get_cycles(;amplitude=amplitude_percent/100, num_steps_per_cycle=Int(round(500*0.5)))
        plot_results!(ax, m, ϵ11; label="ϵ₁₁ = ±$amplitude_percent %")
    end
    CM.axislegend(ax; position=:rb)
    return fig, ax
end

fig1, ax1 = create_example_plot(m_chaboche, [0.15, 0.9])
fig2, ax2 = create_example_plot(m_chaboche, [0.15])

t = "Figure 4. Cyclic response of Chaboche model to three triangular cycles with 𝜖11 = ±0.15% and 𝜖11 = ±0.9%."
CM.Label(fig1[2,1], t)

t2 = "Figure 5. Cyclic response of Chaboche model to three triangular cycles with 𝜖11 = ±0.15% "
CM.Label(fig2[2,1], t2)

CM.save("figure4.png", fig1)
CM.save("figure5.png", fig2)
fig
fig2

end