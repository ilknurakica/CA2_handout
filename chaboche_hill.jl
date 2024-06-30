struct ChabocheHill{T} <: AbstractMaterial
    G::T    # Shear modulus
    K::T    # Bulk modulus
    Y::T   # initial Yield limit
    Hiso::T 
    κ∞::T
    Hkin::T
    β∞::T
    G_Hill::T
    H_Hill::T
    F_Hill::T
end

"""
    Chaboche(; G, K, Y, Hiso, κ∞, Hkin, β∞)

Keyword constructor for the Chaboche struct
"""
ChabocheHill(; G, K, Y, Hiso, κ∞, Hkin, β∞, G_Hill, H_Hill, F_Hill) = ChabocheHill(G, K, Y, Hiso, κ∞, Hkin, β∞, G_Hill, H_Hill, F_Hill)

# Define state varibale struct
struct ChabocheHillState{T} <: AbstractMaterialState
    ϵp::SymmetricTensor{2, 3, T, 6}  # 6 is the number of independent components in the tensor
    β::SymmetricTensor{2, 3, T, 6}
    κ::T
    
    end

MaterialModelsBase.initial_material_state(::Chaboche) = ChabocheHillState(zero(SymmetricTensor{2,3}), zero(SymmetricTensor{2,3}), 0.0)

# Elastic functions
function elastic_stress(m::ChabocheHill, ϵe::SymmetricTensor{2,3})
    return 2*m.G*dev(ϵe) + 3*m.K*vol(ϵe)
end

function elastic_stiffness(m::ChabocheHill)
    I2 = one(SymmetricTensor{2,3})
    I4sym = symmetric(otimesu(I2,I2))
    I4vol = I2⊗I2/3
    return 2*m.G*(I4sym-I4vol) + 3*m.K*I4vol
end

function extract_unknowns(::ChabocheHill, x)
    ϵp = frommandel(SymmetricTensor{2, 3}, x[1:6])  
    β = frommandel(SymmetricTensor{2, 3}, x[7:12]) 
    return ϵp, β, x[13], x[14]
end

function sigma_from_x(m::ChabocheHill, ϵ, x)
    ϵp, _, _, _ = extract_unknowns(m, x)
    return elastic_stress(m, ϵ-ϵp)
end

function hill_eff_stress(σ, m)
    # Create the matrix C


    return (σ ⋅ C ⋅ σ)^0.5
    

end

function Hill48(m::ChabocheHill, σ_1::SymmetricTensor{2,3}, σ_2::SymmetricTensor{2,3})
    σ_eff = 


    function vonmises(σ::SymmetricTensor{2,3})
        σ_dev = dev(σ)
        return sqrt((3.0/2.0) * (σ_dev ⊡ σ_dev))
    end