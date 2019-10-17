# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

dependencies = [
    "https://github.com/TeroFrondelius/umat_binaries_builder/releases/download/v0.5.0/build_umat_binaries.v0.5.0.jl"
]

for build_script in dependencies
    script_name = split(build_script,"/")[end]
    include(download(build_script,script_name))
end
