#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

function test_poisson(x, y, p, bc)
	bounds = (p.o, p.n)

	A    = x -> divergence(grad(x))
	B(λ) = x -> A(x) - λ*x
	
	# A f = b
	f = Field(p.n, div_stags, bc)
	b = Field(p.n, div_stags)
	
	# helpers
	h = Field(p.n, div_stags, bc)
	r = Field(p.n, div_stags, bc)

	# modes
	m_1 = Field(p.n, div_stags, bc)
	m_n = Field(p.n, div_stags, bc)
	
	assign!(b, 1, bounds)
	
	# random initial guesses:
	assign!(f,   fieldgen((_...) -> rand()), bounds)
	assign!(m_1, fieldgen((_...) -> rand()), bounds)
	assign!(m_n, fieldgen((_...) -> rand()), bounds)
	
	λ_n =       powerit!(A,      m_n, r; bounds = bounds, nit = 10000, atol = 1e-7)
	λ_1 = λ_n + powerit!(B(λ_n), m_1, r; bounds = bounds, nit = 10000, atol = 1e-7)

	Meta.@show λ_1, λ_n

	ε = chebyshev!(A, f, b; λ = (λ_1, λ_n), v = h, r = r, bounds = bounds, atol = 1e-7)

	display(plot(log10.(ε)))
	readline()

	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
end

function test_elasticity(x, y, p, bc)
	
end

function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	o    =  n .- n .+ 1              # logical origin
	
	# collect all variables local to this function:
	vars = Base.@locals
	
	# and return as a named tuple
	return NamedTuple{Tuple(keys(vars))}(values(vars))
end

function main()
	s = ArgParseSettings()
	
	@add_arg_table s begin
    	"--nb"
			help = "comma-separated list of the number of blocks in each dimension"
			required = true
			arg_type = (NTuple{N, Int} where N)
    end
	
	p = parameters(; parse_args(s)...); pprintln(p)
	
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	#    test_poisson(x, y, p, (Essential(:-, :y, 0), Essential(:-, :x, 0), Essential(:+, :x, 0), Essential(:+, :y, 1)))
	   test_poisson(x, y, p, Essential())
	test_elasticity(x, y, p, ImpermeableFreeSlip())
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end