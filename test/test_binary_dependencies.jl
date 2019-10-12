# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

pkg_dir = dirname(Base.find_package("UMAT"))
usr_dir = joinpath(pkg_dir,"..","deps","usr")
if Sys.iswindows()
    @test isfile(joinpath(usr_dir,"bin","libmises_umat.dll"))
else
    @test isfile(joinpath(usr_dir,"lib","libmises_umat.so"))
end
