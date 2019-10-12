# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT, Test, Materials, Libdl

pkg_dir = dirname(Base.find_package("UMAT"))
if Sys.iswindows()
    lib_dir = joinpath(pkg_dir,"..","deps","usr","bin")
else
    lib_dir = joinpath(pkg_dir,"..","deps","usr","lib")
end

E = 160e9; NU = 0.3; SYIELD = 100
parameters = UmatParameterState([E, NU, SYIELD])
material = UmatMaterial(parameters=parameters,
                        lib_path=joinpath(lib_dir,"libmises_umat." * dlext),
                        behaviour=:mises_umat_)
