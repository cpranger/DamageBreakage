# const Float = Float64

const JULIA_CUDA_USE_BINARYBUILDER = false
const USE_GPU = false
const GPU_ID  = 1

const BLOCK_SIZE = 1 # size of an array block

<|(f, x) = f(x)

verbosity::Int = 999
algo_depth::Int = 0
out::Union{IOStream, Base.DevNull} = devnull

macro algo_step(expr)
	return esc <| quote
		Main.algo_depth += 1
		$expr
		Main.algo_depth -= 1
	end
end

macro verbo_println(expr)
	return esc <| quote
		if Main.algo_depth <= Main.verbosity
			println( ' '^(2*Main.algo_depth), $expr)
		end
		println(out, ' '^(2*Main.algo_depth), $expr)
	end
end

macro verbo_error(expr)
	return esc <| quote
		flush(out)
		error( ' '^max(2*Main.algo_depth-7, 0), $expr)
	end
end

using PrettyPrinting
using Profile
using ArgParse
using Plots
using Printf

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

includet("../algorithm/visu.jl")
includet("../algorithm/lanczos.jl")
includet("../algorithm/powerit.jl")
includet("../algorithm/rayleighquotientit.jl")
includet("../algorithm/newtonit.jl")
includet("../algorithm/chebyshev.jl")
includet("../algorithm/cg.jl")
includet("../algorithm/gershgorin.jl")
includet("../algorithm/tr-bdf2.jl")