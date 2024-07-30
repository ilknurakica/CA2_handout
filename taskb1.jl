import CairoMakie as CM
function main_yield_surface()
function yield_surface(F,G,H,M,σ, α)
σ_1 = σ*(cos(α))^2
σ_2 = σ*(sin(α))^2
σ_12 = σ*(cos(α))*sin(α)
    return sqrt(F*(σ_1 - σ_2)^2 + G*σ_1^2 + H*σ_2^2 + M*σ_12^2) - 1

end 

# define array of angles
angles = collect(range(0, stop=π/2, length=30))

# define material parameters

G = 0.625
H = 0.125
M = 1.725*2
F = 0.375

F_vm=0.5
G_vm=0.5
H_vm=0.5
M_vm=3.0

# Define the applied stress
σ = 1

# Initialize an array to store the computed effective stress values
effective_stress = Float64[]
eff_stress_vm = Float64[]

# Compute the effective stress for each angle
for α in angles
    push!(effective_stress, yield_surface(F, G, H, M, σ, α))
    push!(eff_stress_vm, yield_surface(F_vm, G_vm, H_vm, M_vm, σ, α))
end

# Initialize a figure
fig = CM.Figure()

# Create an axis for plotting
ax = CM.Axis(fig[1, 1], title="Figure 7. Yield surface for angles between 0 and pi/2 in case of uniaxial in plane loading of a sheet.
", xlabel="Angle (radians)", ylabel="Yield Surface")

# Plot the computed effective stress values
CM.lines!(ax, angles, effective_stress, label="Anisotropic")
CM.lines!(ax, angles, eff_stress_vm, label="Isotropic")

# Add a legend
CM.axislegend(ax; position=:rb)
CM.save("figure7.png", fig)
# Display the plot
CM.display(fig)

end