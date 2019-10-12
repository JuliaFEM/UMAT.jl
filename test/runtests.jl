# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

using UMAT
using Test
using Suppressor


@testset "UMAT.jl" begin
    @testset "Binary dependencies" begin
        include("test_binary_dependencies.jl")
    end
    @testset "libelastic.so" begin
        include("test_elastic.jl")
    end
    @testset "UmatMaterial" begin
        include("test_umatmaterial.jl")
    end
end
