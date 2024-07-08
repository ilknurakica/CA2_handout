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


G_hill = 0.125
H_hill = 0.625
F_hill = 0.375
K_hill = 1.500
L_hill = 1.500
M_hill = 1.725

C = [
  F + G  -F    -G    0    0    0
  -F    F + H  -H    0    0    0
  -G    -H    G + H  0    0    0
  0     0     0     K/2  0    0
  0     0     0     0    L/2  0
  0     0     0     0    0    M/2
]

C2 = [
  F_hill + G_hill  -F_hill    -G_hill    0    0    0
  -F_hill    F_hill + H_hill  -H_hill    0    0    0
  -G_hill    -H_hill    G_hill + H_hill  0    0    0
  0     0     0     K_hill/2  0    0
  0     0     0     0    L_hill/2  0
  0     0     0     0    0    M_hill/2
]
C = frommandel(SymmetricTensor{4,3}, C)
C2 = frommandel(SymmetricTensor{4,3}, C2)
a = rand(SymmetricTensor{2, 3})
b_c = sqrt(a ⊡ C ⊡ a)
b = sqrt(a ⊡ C2 ⊡ a)
isapprox(b_c, b)
print(b_c)
print(b)




