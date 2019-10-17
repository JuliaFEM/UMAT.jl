# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT, Test, Materials, Libdl

pkg_dir = dirname(Base.find_package("UMAT"))
if Sys.iswindows()
    lib_dir = joinpath(pkg_dir,"..","deps","usr","bin")
else
    lib_dir = joinpath(pkg_dir,"..","deps","usr","lib")
end

material = UmatMaterial(NTENS=4, NSTATV=0, NPROPS=0,
                        lib_path=joinpath(lib_dir,"libelastic." * dlext))

strains = [1.0e-3, 0.0, 0.0, 0.0]
ref_stress = [282.6923171355886, 121.1538570784259, 121.1538570784259, 0.0]
material.ddrivers = UmatDriverState(NTENS=4, STRAN = strains)
UMAT.integrate_material!(material)
update_material!(material)

for i in 1:4
    @test isapprox(material.variables.STRESS[i],ref_stress[i],atol=sqrt(eps()))
end
