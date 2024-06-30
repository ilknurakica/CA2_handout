import CairoMakie as CM

function yield_surface(N,G,H,σ, α)
σ_1 = σ*cos(α)
σ_12 = σ*sin(α)

    return sqrt((G+H)*σ_1^2 + 2*N*σ_12^2) - 1

end 

# define array of angles
angles = collect(range(0, stop=π/2, length=30))

# define material parameters

G = 0.625
H = 0.375
N = 1.725

# Define the applied stress
σ = 1

# Initialize an array to store the computed effective stress values
effective_stress = Float64[]

# Compute the effective stress for each angle
for α in angles
    push!(effective_stress, yield_surface(N, G, H, σ, α))
end

# Initialize a figure
fig = CM.Figure()

# Create an axis for plotting
ax = CM.Axis(fig[1, 1], title="Yiled surface for Angles between 0 and π/2", xlabel="Angle (radians)", ylabel="Yield Surface")

# Plot the computed effective stress values
CM.lines!(ax, angles, effective_stress, label="Yield Surface")

# Add a legend
#CM.legend!(ax)

# Display the plot
CM.display(fig)