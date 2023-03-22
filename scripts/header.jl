# const Float = Float64

const JULIA_CUDA_USE_BINARYBUILDER = false
const USE_GPU = false
const GPU_ID  = 1

const BLOCK_SIZE = 1 # size of an array block

using Revise

<|(f, x) = f(x)

using PrettyPrinting
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

includet("../algorithm/algorithm.jl")
includet("../algorithm/visu.jl")

using .Algo