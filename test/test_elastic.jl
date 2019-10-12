# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT, Test, Materials, Libdl

pkg_dir = dirname(Base.find_package("UMAT"))
if Sys.iswindows()
    lib_dir = joinpath(pkg_dir,"..","deps","usr","bin")
else
    lib_dir = joinpath(pkg_dir,"..","deps","usr","lib")
end

material = UmatMaterial(lib_path=joinpath(lib_dir,"libelastic." * dlext),
                        behaviour=:umat_)

strains = [1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0]
material.ddrivers = UmatDriverState(STRAN = strains)
UMAT.integrate_material!(material)
update_material!(material)


@test isapprox(material.variables.STRESS[1],210.)
for i in 2:6
    @test isapprox(material.variables.STRESS[i],0.0)
end
