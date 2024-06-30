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


