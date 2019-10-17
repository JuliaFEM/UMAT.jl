# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

cmake_minimum_required(VERSION 3.5)
project(UMAT.jl/umat_models)

set(VERSION 0.5.0)
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-m64 -fdefault-real-8 -fdefault-integer-8 -free")

add_library(gurson_porous_plasticity SHARED gurson_porous_plasticity.f90)
target_link_libraries(gurson_porous_plasticity openblas64_)

add_library(drucker_prager_plasticity SHARED drucker_prager_plasticity.f90 drucker_prager_plasticity_deriv.f90)


install(TARGETS gurson_porous_plasticity DESTINATION lib)
install(TARGETS drucker_prager_plasticity DESTINATION lib)
