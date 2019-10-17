# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

material = DruckerPragerMaterial()

# Define strain history
e11 = 0.5.*vcat(Array(range(0, stop=1.5e-2, length=30)), Array(range(1.5e-2, stop=-1.5e-2, length=30)), Array(range(-1.5e-2, stop=3e-2, length=40)))
e22 = -e11/2
e33 = e22
strains = [[e11[i], e22[i], e33[i], 0.0, 0.0, 0.0] for i in 1:100]

s11 = [material.variables.STRESS[1]]
s22 = [material.variables.STRESS[2]]
s33 = [material.variables.STRESS[3]]

dtime = 0.01

for i=2:100
    dstrain = strains[i]-strains[i-1]
    material.ddrivers = UmatDriverState(NTENS=6, STRAN=dstrain)
    UMAT.integrate_material!(material)
    push!(s11, material.variables.STRESS[1])
    push!(s22, material.variables.STRESS[2])
    push!(s33, material.variables.STRESS[3])
    Materials.update_material!(material)
end
