# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT, Test

pkg_dir = dirname(Base.find_package("UMAT"))
usr_dir = joinpath(pkg_dir,"..","deps","usr")

if !Sys.iswindows()
    material = UmatMaterial(lib_path=joinpath(usr_dir,"lib","libmises_umat.so"),
                            behaviour=:mises_umat_)
end
