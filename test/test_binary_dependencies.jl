# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using Test, Libdl

pkg_dir = dirname(Base.find_package("UMAT"))
if Sys.iswindows()
    lib_dir = joinpath(pkg_dir,"..","deps","usr","bin")
else
    lib_dir = joinpath(pkg_dir,"..","deps","usr","lib")
end

@test isfile(joinpath(lib_dir,"libelastic." * dlext))
@test isfile(joinpath(lib_dir,"libisotropic_plast_exp." * dlext))
@test isfile(joinpath(lib_dir,"libisotropic_plast_imp." * dlext))
@test isfile(joinpath(lib_dir,"libmises_umat." * dlext))
