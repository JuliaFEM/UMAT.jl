# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

module UMAT

using Libdl
using Materials
using Parameters
using LinearAlgebra
using Tensors

pkg_dir = dirname(Base.find_package("UMAT"))
if Sys.iswindows()
    lib_dir = joinpath(pkg_dir,"..","deps","usr","bin")
else
    lib_dir = joinpath(pkg_dir,"..","deps","usr","lib")
end


"""
Variables updated by UMAT routine.
"""
@with_kw mutable struct UmatVariableState <: AbstractMaterialState
    NTENS :: Int
    NSTATV :: Int = zero(Int)
    DDSDDE :: Array{Float64,2} = zeros(Float64, NTENS, NTENS)
    STRESS :: Array{Float64,1} = zeros(Float64, NTENS)
    STATEV :: Array{Float64,1} = zeros(Float64, NSTATV)
    SSE :: Array{Float64,1} = zeros(Float64, 1)
    SPD :: Array{Float64,1} = zeros(Float64, 1)
    SCD :: Array{Float64,1} = zeros(Float64, 1)
    RPL :: Array{Float64,1} = zeros(Float64, 1)
    DDSDDT :: Array{Float64,1} = zeros(Float64, NTENS)
    DRPLDE :: Array{Float64,1} = zeros(Float64, NTENS)
    DRPLDT :: Array{Float64,1} = zeros(Float64, 1)
    PNEWDT :: Array{Float64,1} = ones(Float64, 1)
    # TODO Check that shear stress components are in correct order
    # https://github.com/KristofferC/Tensors.jl/blob/e3a67612f38124c15d77ea0a2a21d0175d0a32d8/src/voigt.jl#L1
    stress = fromvoigt(SymmetricTensor{2,convert(Int,NTENS/2)}, STRESS)
    jacobian = fromvoigt(SymmetricTensor{4,convert(Int,NTENS/2)}, DDSDDE)
end

"""
Material parameters in order that is specific to chosen UMAT.
"""
@with_kw mutable struct UmatParameterState <: AbstractMaterialState
    NPROPS :: Int = zero(Int)
    PROPS :: Array{Float64,1} = zeros(Float64, NPROPS)
end

"""
Variables passed in for information.
These drive evolution of the material state.
"""
@with_kw mutable struct UmatDriverState <: AbstractMaterialState
    NTENS :: Int
    STRAN :: Array{Float64,1} = zeros(Float64, NTENS)
    TIME :: Array{Float64,1} = zeros(Float64, 2)
    TEMP :: Float64 = zero(Float64)
    PREDEF :: Float64 = zero(Float64)
    # TODO Check that shear strain components are in correct order
    # https://github.com/KristofferC/Tensors.jl/blob/e3a67612f38124c15d77ea0a2a21d0175d0a32d8/src/voigt.jl#L1
    strain = fromvoigt(SymmetricTensor{2,convert(Int,NTENS/2)}, STRAN; offdiagscale = 2.0)
    time = TIME[2]
end

"""
Other Abaqus UMAT variables passed in for information.
"""
@with_kw mutable struct UmatOtherState <: AbstractMaterialState
    CMNAME :: String = "U" # Ptr{Cuchar} # = 'U'
    NDI :: Int = 3
    NSHR :: Int = 3
    NTENS :: Int  # NDI + NSHR
    NSTATV :: Int
    NPROPS :: Int
    COORDS :: Array{Float64,1} = zeros(Float64, 3)
    DROT :: Array{Float64,2} = I + zeros(Float64, 3, 3)
    CELENT :: Float64 = zero(Float64)
    DFGRD0 :: Array{Float64,2} = I + zeros(Float64, 3, 3)
    DFGRD1 :: Array{Float64,2} = I + zeros(Float64, 3, 3)
    NOEL :: Int = zero(Int)
    NPT :: Int = zero(Int)
    LAYER :: Int = zero(Int)
    KSPT :: Int = zero(Int)
    JSTEP :: Array{Float64,1} = zeros(Int, 3)
    KSTEP :: Int = zero(Int)
    KINC :: Int = zero(Int)
end

"""
UMAT material structure.

`lib_path` is the path to the compiled shared library.
`behaviour` is the function name to call from the shared library
    generally this is compiler dependent, by default `umat_`.
    MFront Abaqus interface produces specific name, e.g. `ELASTICITY_3D`.
"""
@with_kw mutable struct UmatMaterial <: AbstractMaterial
    NTENS :: Int
    NSTATV :: Int = zero(Int)
    NPROPS :: Int = zero(Int)

    drivers :: UmatDriverState = UmatDriverState(NTENS=NTENS)
    ddrivers :: UmatDriverState = UmatDriverState(NTENS=NTENS)
    variables :: UmatVariableState = UmatVariableState(NTENS=NTENS, NSTATV=NSTATV)
    variables_new :: UmatVariableState = UmatVariableState(NTENS=NTENS, NSTATV=NSTATV)
    parameters :: UmatParameterState = UmatParameterState(NPROPS=NPROPS)
    dparameters :: UmatParameterState = UmatParameterState(NPROPS=NPROPS)

    umat_other :: UmatOtherState = UmatOtherState(NTENS=NTENS, NSTATV=NSTATV, NPROPS=NPROPS)

    lib_path :: String
    behaviour :: Symbol = :umat_

end


"""
Wrapper function to ccall the UMAT subroutine from the compiled shared library.
"""
function call_umat!(func_umat::Symbol, lib_path::String, STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
        STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,
        NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,
        NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

    lib_umat = Libdl.dlopen(lib_path) # Open the library explicitly.
    sym_umat = Libdl.dlsym(lib_umat, func_umat) # Get a symbol for the umat function to call.

    ccall(sym_umat, Nothing,
        (Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
        Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ptr{Cuchar},#Csize_t,
        Ref{Int},Ref{Int},
        Ref{Int},Ref{Int},Ref{Float64},Ref{Int},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
        Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int}),
        STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
        STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,#sizeof(CMNAME),
        NDI,NSHR,
        NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,
        NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

    Libdl.dlclose(lib_umat) # Close the library explicitly.

    return Nothing
end

"""
Calls UMAT and updates MaterialVariableState writing the result to material.variables_new
"""
function Materials.integrate_material!(material::UmatMaterial)
    @unpack CMNAME,NDI,NSHR,NTENS,NSTATV,NPROPS,COORDS,DROT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC = material.umat_other
    @unpack PROPS = material.parameters
    STRAN, TIME, TEMP, PREDEF = material.drivers.STRAN, material.drivers.TIME, material.drivers.TEMP, material.drivers.PREDEF
    DSTRAN, DTIME, DTEMP, DPRED = material.ddrivers.STRAN, material.drivers.TIME[1], material.ddrivers.TEMP, material.ddrivers.PREDEF
    @unpack STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,PNEWDT = material.variables

    call_umat!(material.behaviour, material.lib_path, STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT,
        STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR,
        NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,
        NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    variables_new = UmatVariableState(NSTATV=NSTATV,NTENS=NTENS,STRESS=STRESS,STATEV=STATEV,DDSDDE=DDSDDE,SSE=SSE,SPD=SPD,SCD=SCD,RPL=RPL,DDSDDT=DDSDDT,DRPLDE=DRPLDE,DRPLDT=DRPLDT,PNEWDT=PNEWDT)
    material.variables_new = variables_new
end

function Materials.reset_material!(material::UmatMaterial)
    material.ddrivers = UmatDriverState(NTENS=material.NTENS)
    material.dparameters = UmatParameterState(NPROPS=material.NPROPS)
    material.variables_new = UmatVariableState(NTENS=material.NTENS, NSTATV=material.NSTATV)
    return nothing
end

include("gurson_model.jl")
include("druckerprager_model.jl")

export UmatMaterial, UmatDriverState, UmatParameterState, UmatVariableState
export UmatOtherState, GursonMaterial, DruckerPragerMaterial, integrate_material!

# Re-export update_material! and uniaxial_increment! from Materials.jl for convenience
export update_material!, uniaxial_increment!
end # module
