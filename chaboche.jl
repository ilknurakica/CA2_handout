
# Define the Chaboche struct
struct Chaboche{T} <: AbstractMaterial
    G::T    # Shear modulus
    K::T    # Bulk modulus
    Y::T   # initial Yield limit
    Hiso::T 
    κ∞::T
    Hkin::T
    β∞::T
end

"""
    Chaboche(; G, K, Y, Hiso, κ∞, Hkin, β∞)

Keyword constructor for the Chaboche struct
"""
Chaboche(; G, K, Y, Hiso, κ∞, Hkin, β∞) = Chaboche(G, K, Y, Hiso, κ∞, Hkin, β∞)

# Define state varibale struct
struct ChabocheState{T} <: AbstractMaterialState
    ϵp::SymmetricTensor{2, 3, T, 6}  # 6 is the number of independent components in the tensor
    β::SymmetricTensor{2, 3, T, 6}
    κ::T
    end

MaterialModelsBase.initial_material_state(::Chaboche) = ChabocheState(zero(SymmetricTensor{2,3}), zero(SymmetricTensor{2,3}), 0.0)

# Elastic functions
function elastic_stress(m::Chaboche, ϵe::SymmetricTensor{2,3})
    return 2*m.G*dev(ϵe) + 3*m.K*vol(ϵe)
end

function elastic_stiffness(m::Chaboche)
    I2 = one(SymmetricTensor{2,3})
    I4sym = symmetric(otimesu(I2,I2))
    I4vol = I2⊗I2/3
    return 2*m.G*(I4sym-I4vol) + 3*m.K*I4vol
end

function extract_unknowns(::Chaboche, x)
    ϵp = frommandel(SymmetricTensor{2, 3}, x[1:6])  
    β = frommandel(SymmetricTensor{2, 3}, x[7:12]) 
    return ϵp, β, x[13], x[14]
end

function sigma_from_x(m::Chaboche, ϵ, x)
    ϵp, _, _, _ = extract_unknowns(m, x)
    return elastic_stress(m, ϵ-ϵp)
end

# main function
function MaterialModelsBase.material_response(
    m::Chaboche, ϵ::SymmetricTensor{2,3}, 
    old_state::ChabocheState, args...)
    σ_trial = elastic_stress(m, ϵ-old_state.ϵp)
    ϕ = sqrt(3/2)*norm(dev(σ_trial) - dev(old_state.β)) - (m.Y + old_state.κ)
    if ϕ < 0
        # Elastic
        return σ_trial, elastic_stiffness(m), old_state
    else
        # Plastic
        rf!(r_,x_) = residual!(r_, x_, m, ϵ, old_state)
        x0 = zeros(14)               # Allocate unknowns
        # Set starting guess
        tomandel!(x0[1:6], old_state.ϵp)  
        tomandel!(x0[7:12], old_state.β)
        x0[14] = old_state.κ
        x, ∂r∂x, converged = newtonsolve(rf!, x0)
        if converged
            σ = sigma_from_x(m, ϵ, x)
            dσdϵ = calculate_ats(m, x, ϵ, old_state, ∂r∂x) #???
            ϵp, β, _, κ = extract_unknowns(m, x)
            new_state = ChabocheState(ϵp, β, κ) 
            return σ, dσdϵ, new_state
        else
            throw(ErrorException("PerfectPlasticity did not converge"))
        end
    end
end


function calculate_ats(m::Chaboche, x::Vector, ϵ::SymmetricTensor{2,3}, old_state, ∂r∂x::Matrix)
    # dσdϵ = ∂σ∂ϵ + ∂σ∂x ∂x∂ϵ
    # Problem: x is an implicit function of ϵ
    # drdϵ = 0 = ∂r∂ϵ + ∂r∂x ∂x∂ϵ => ∂x∂ϵ = -∂r∂x\∂r∂ϵ

    ∂σ∂ϵ = elastic_stiffness(m)
    rf_ϵ!(r_, ϵv_) = residual!(r_, x, m, frommandel(SymmetricTensor{2,3}, ϵv_), old_state)
    ∂r∂ϵ = ForwardDiff.jacobian(rf_ϵ!, zeros(14), tomandel(ϵ))
    ∂x∂ϵ = -∂r∂x\∂r∂ϵ
    
    ∂σ∂x = ForwardDiff.jacobian(x_ -> tomandel(sigma_from_x(m, ϵ, x_)), x)
    return ∂σ∂ϵ + frommandel(SymmetricTensor{4,3}, ∂σ∂x*∂x∂ϵ)
end





function residual!(r::Vector, x::Vector, m::Chaboche, ϵ::SymmetricTensor{2,3}, old_state::ChabocheState)
    ϵp, β, Δλ, κ = extract_unknowns(m, x)
    σ = elastic_stress(m, ϵ-ϵp)
    ν = (3/2)*(dev(σ) - dev(β))/vonmises(σ - β)
    #ν = gradient(vonmises, σ)
    
    r1 = (ϵp - old_state.ϵp) - Δλ*ν
    # k = gradient(r1,ϵp )
    # k= k
    r2 = β - old_state.β - Δλ*(2/3)*m.Hkin*(ν - (3/2)*(β/m.β∞))
    r3 = vonmises(σ - β) - (m.Y + κ)
    r4 = κ - old_state.κ - Δλ*m.Hiso*(1 - κ/m.κ∞)
    
    
    r[1:6] .= tomandel(r1)   # Store the Mandel representation of r1 in the first 6 elements
    r[7:12] .= tomandel(r2)
    r[13] = r3
    r[14] = r4
end