# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT, Test, Materials, Libdl, Tensors
const abaqus = [1 4 6; 0 2 5; 0 0 3]

pkg_dir = dirname(Base.find_package("UMAT"))
if Sys.iswindows()
    lib_dir = joinpath(pkg_dir,"..","deps","usr","bin")
else
    lib_dir = joinpath(pkg_dir,"..","deps","usr","lib")
end

material = UmatMaterial(NTENS=4, NSTATV=0, NPROPS=0,
                        lib_path=joinpath(lib_dir,"libelastic." * dlext))

strains = fromvoigt(SymmetricTensor{2,2},[1.0e-3, 0.0, 0.0, 0.0]; offdiagscale = 2.0, order=abaqus)
ref_stress = [282.6923171355886, 121.1538570784259, 121.1538570784259, 0.0]
material.ddrivers = UmatDriverState(NTENS=4, strain=strains)
UMAT.integrate_material!(material)
update_material!(material)


for i in 1:3
    @test_broken isapprox(tovoigt(material.variables.stress, order=abaqus)[i],ref_stress[i],atol=sqrt(eps()))
end
