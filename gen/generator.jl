#! /bin/bash julia --project generator.jl
using Pkg
using Pkg.Artifacts
using Clang.Generators
using Clang.Generators.JLLEnvs
using SuiteSparse_jll
using MPFR_jll
using GMP_jll

cd(@__DIR__)

# headers
stdlib_toml = joinpath(dirname(pathof(SuiteSparse_jll)), "..", "StdlibArtifacts.toml")
SuiteSparse_dir = Pkg.Artifacts.ensure_artifact_installed("SuiteSparse", stdlib_toml)

ss_dir = joinpath(SuiteSparse_dir, "include") |> normpath

util_h = joinpath(@__DIR__, "SPEX_Util.h")
leftlu_h = joinpath(@__DIR__, "SPEX_Left_LU.h")
@assert isfile(util_h)
@assert isfile(leftlu_h)

# load common option
options = load_options(joinpath(@__DIR__, "generator.toml"))

# run generator for all platforms
args = get_default_args()
push!(args, "-I$ss_dir")
push!(args, "-I~/julia/usr/include")
header_files = [util_h, leftlu_h]
ctx = create_context(header_files, args, options)
build!(ctx)