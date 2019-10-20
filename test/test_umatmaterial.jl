# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT, Test, Materials, Libdl, Tensors

pkg_dir = dirname(Base.find_package("UMAT"))
if Sys.iswindows()
    lib_dir = joinpath(pkg_dir,"..","deps","usr","bin")
else
    lib_dir = joinpath(pkg_dir,"..","deps","usr","lib")
end

# E = 160e9; NU = 0.3; SYIELD = 100
PROPS = [160.e9, 0.3, 100.]
material = UmatMaterial(NTENS=6, NSTATV=0, NPROPS=3,
                        parameters=UmatParameterState(NPROPS=3,PROPS=PROPS),
                        lib_path=joinpath(lib_dir,"libmises_umat." * dlext),
                        behaviour=:mises_umat_)

dtime = 0.1
stresses = [copy(tovoigt(material.variables.stress))]
dstrain11 = 1e-3*dtime
dtimes = [dtime, dtime, dtime, dtime, 1.0]
dstrains11 = [dstrain11, dstrain11, dstrain11, -dstrain11, -4*dstrain11]

expected_S11 = [50. 100. 150. 100. -100.]
expected_stress = zeros(6)

expected_strains = [[0.00025, -7.5e-5, -7.5e-5, 0.0, 0.0, 0.0],
                    [0.0005, -0.00015, -0.00015, 0.0, 0.0, 0.0],
                    [0.00075, -0.000225, -0.000225, 0.0, 0.0, 0.0],
                    [0.0005, -0.00015, -0.00015, 0.0, 0.0, 0.0],
                    [-0.0005, 0.00015, 0.00015, 0.0, 0.0, 0.0]]


for i in 1:length(dtimes)
    dstrain11 = dstrains11[i]
    dtime = dtimes[i]
    uniaxial_increment!(material, dstrain11, dtime)
    update_material!(material)
    expected_stress[1] = expected_S11[i]
    @info(tovoigt(mat.variables.stress), stresses_expected[i])
    @info(tovoigt(material.drivers.strain; offdiagscale=2.0))
    #@test isapprox(tovoigt(material.variables.stress), expected_stress)
    #@test isapprox(tovoigt(material.drivers.strain; offdiagscale=2.0), expected_strains[i])
end
