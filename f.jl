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
  m.F + m.G  -m.F    -m.G    0    0    0
  -m.F    m.F + m.H  -m.H    0    0    0
  -m.G    -m.H    m.G + m.H  0    0    0
  0     0     0     m.K/2  0    0
  0     0     0     0    m.L/2  0
  0     0     0     0    0    m.M/2
]

C = frommandel(SymmetricTensor{4,3}, C)
a = rand(SymmetricTensor{2, 3})
b = vonmises(a)


$ git config --global user.email ilknurakica@gmail.com