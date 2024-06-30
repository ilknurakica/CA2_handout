using Tensors
#using MaterialModelsBase
# using Newton


function vonmises(σ::SymmetricTensor{2,3})
    σ_dev = dev(σ)
    return sqrt((3.0/2.0) * (σ_dev ⊡ σ_dev))
end


  
F=0.5
G=0.5
H=0.5
K=3
L=3
M=3
C = [
  F + G  -F    -G    0    0    0
  -F    F + H  -H    0    0    0
  -G    -H    G + H  0    0    0
  0     0     0     K/2  0    0
  0     0     0     0    L/2  0
  0     0     0     0    0    M/2
]

C = frommandel(SymmetricTensor{4,3}, C)
a = rand(SymmetricTensor{2, 3})
b_c = sqrt(a ⊡ C ⊡ a)
b = vonmises(a)

isapprox(b_c, b)


"""
    function uniaxial_stress(m::AbstractMaterial, ϵ11::Vector, t::Vector)    
    
For a time history of the ``\\epsilon_{11}`` component, return 
history of the ``\\sigma_{11}`` component assuming a uniaxial stress
state. Requires that `initial_material_state(m)` and `material_response(m, ...)`
have been defined. 
"""
"""
function uniaxial_stress(m::AbstractMaterial, ϵ11::Vector, t::Vector)
    state = initial_material_state(m)
    σ11 = zeros(length(ϵ11))
    ϵ = zero(SymmetricTensor{2,3})
    for i in 2:length(ϵ11)
        Δt = t[i]-t[i-1]
        ϵ_guess = SymmetricTensor{2,3}((k,l)->k==l==1 ? ϵ11[i] : ϵ[k,l])
        ϵ, σ, state = uniaxial_stress_iterations(m, state, Δt, ϵ_guess)
        σ11[i] = σ[1,1]
    end
    return σ11 
end
"""


