# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT, Test, DelimitedFiles, Tensors

dtime = 0.25
mat = GursonMaterial()
times = [mat.drivers.time]
stresses = [copy(tovoigt(mat.variables.stress))]
dstrain11 = 1e-3*dtime
dtimes = [dtime, dtime, dtime, dtime, 1.0]
dstrains11 = [dstrain11, dstrain11, dstrain11, -dstrain11, -4*dstrain11]
# stresses_expected = [[38.4091, 0.0, 0.0, 0.0, 0.0, 0.0],
#                      [100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
#                      [100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
#                      [50.0, 0.0, 0.0, 0.0, 0.0, 0.0],
#                      [-100.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
# strains_expected = [[dstrain11, -0.3*dstrain11, -0.3*dstrain11, 0.0, 0.0, 0.0],
#                     [2*dstrain11, -0.3*dstrain11*2, -0.3*dstrain11*2, 0.0, 0.0, 0.0],
#                     [3*dstrain11, -0.3*dstrain11*2 - 0.5*dstrain11, -0.3*dstrain11*2 - 0.5*dstrain11, 0.0, 0.0, 0.0],
#                     [2*dstrain11, -0.3*dstrain11 - 0.5*dstrain11, -0.3*dstrain11 - 0.5*dstrain11, 0.0, 0.0, 0.0],
#                     [-2*dstrain11, 0.3*dstrain11*2, 0.3*dstrain11*2, 0.0, 0.0, 0.0]]
#for i in 1:length(dtimes)
    #dstrain11 = dstrains11[i]
    #dtime = dtimes[i]
    uniaxial_increment!(mat, dstrain11, dtime)
    update_material!(mat)
    # push!(times, mat.drivers.time)
    # push!(stresses, copy(tovoigt(mat.variables.stress)))
    #@info(tovoigt(mat.variables.stress), stresses_expected[i])
    #@test isapprox(tovoigt(mat.variables.stress), stresses_expected[i])
    #@test isapprox(tovoigt(mat.drivers.strain; offdiagscale=2.0), strains_expected[i])
#end

# Define strain history
#e11 = 0.5.*vcat(Array(range(0, stop=1.5e-2, length=30)), Array(range(1.5e-2, stop=-1.5e-2, length=30)), Array(range(-1.5e-2, stop=3e-2, length=40)))
#e22 = -e11/2
#e33 = e22
#strains = [[e11[i], e22[i], e33[i], 0.0, 0.0, 0.0] for i in 1:100]

#s11 = [material.variables.STRESS[1]]
#s22 = [material.variables.STRESS[2]]
#s33 = [material.variables.STRESS[3]]

# dtime = 0.01
#
# for i=2:100
#     dstrain = strains[i]-strains[i-1]
#     material.ddrivers = UmatDriverState(NTENS=6, STRAN=dstrain)
#     UMAT.integrate_material!(material)
#     push!(s11, material.variables.STRESS[1])
#     push!(s22, material.variables.STRESS[2])
#     push!(s33, material.variables.STRESS[3])
#     Materials.update_material!(material)
# end
#
# for sig in ["s11", "s22", "s33"]
#     ref_stress = readdlm("test_gurson/ref_stresses_" * sig *".txt")
#     for (comp,ref) in zip(eval(Symbol(sig)), ref_stress)
#         @test isapprox(comp,ref,atol=sqrt(eps()))
#     end
# end
