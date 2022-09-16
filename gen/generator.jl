#! /bin/bash julia --project generator.jl
using Pkg
using Pkg.Artifacts
using Clang.Generators
using SuiteSparse_jll

cd(@__DIR__)

# headers

spex_inc = "/Users/raykimm/SPEX/include/"

util_h = joinpath(spex_inc, "SPEX_Util.h")
leftlu_h = joinpath(spex_inc, "SPEX_Left_LU.h")
@assert isfile(util_h)
@assert isfile(leftlu_h)

# load common option
options = load_options(joinpath(@__DIR__, "generator.toml"))

# run generator for all platforms
args = get_default_args()
push!(args, "-I$spex_inc")
push!(args, "-isystem/opt/homebrew/include")
header_files = [util_h, leftlu_h]
ctx = create_context(header_files, args, options)
build!(ctx)