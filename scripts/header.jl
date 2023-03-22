const Float = Float64

const JULIA_CUDA_USE_BINARYBUILDER = false
const USE_GPU = false
const GPU_ID  = 1

const BLOCK_SIZE = 1 # size of an array block

using Revise

using ParallelStencil
# using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
	@init_parallel_stencil(CUDA, Float, 2)
	CUDA.device!(GPU_ID)
else
	@init_parallel_stencil(Threads, Float, 2)
end
using JLD
using Printf
using Adapt
using PrettyPrinting
using Plots
using Profile

using ArgParse

function ArgParse.parse_item(::Type{(NTuple{N, Int} where N)}, x::AbstractString)
	return Tuple <| map(ss -> parse(Int, ss), split(x, ','))
end

function parse_args(s)
    args = ArgParse.parse_args(ARGS, s, as_symbols = true)
	vars = Tuple <| keys(args)
	vals = Tuple <| values(args)
	return (; zip(vars, vals)...)
end

using StaggeredKernels

include("../algorithm/algorithm.jl")
include("../algorithm/visu.jl")

using .Algo