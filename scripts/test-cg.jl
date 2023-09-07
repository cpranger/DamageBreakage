#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

gen_rand() = fieldgen((_...) -> rand())

(gen_rand(v::Tensor{S, NamedTuple{N}}) where {S, N}) = (; zip(N, gen_rand() for n in N...)...)

function test_poisson(p)
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	bounds = (p.o, p.n)

	f   = Field(p.n, div_stags)
	h_1 = deepcopy(f)
	h_2 = deepcopy(f)
	r   = deepcopy(f)
	b   = deepcopy(f)
	
	mode = Mode(f, Essential())
	println("λ* = ($(-mode[-1,-1].val^2), $(-mode[1,1].val^2))")
	
	R = x -> (
		divergence(grad(x)) - 1/1000,
		(
			Essential(:-, :y, FD(x, :y)),
			Essential(:-, :x, -x),
			Essential(:+, :x, -x),
			Essential(:+, :y, -x + 1)
		)
	)

	A = linearize(R, 0)
	assign!(b, -R(0), bounds)

	assign!(f, (0, R(0)[2]), bounds)
	(λ, Λ, ε) = cg!(A, f, b; r = r, p = h_1, q = h_2, bounds = bounds, rtol = 1e-8, λtol = 1e-2, minit = 100, maxit = 1000)
	display <| plot(log10.(ε))
	readline()
	
	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	display <| plot(plt1, plt2; layout = (1, 2))
end

function test_elasticity(p)
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	bounds = (p.o, p.n)

	A    = x -> divergence(symgrad(x))

	# A f = b
	f   = Vector(p.n, motion_stags)
	b   = Vector(p.n, motion_stags)
	
	bc_expl_mode = (
		x = (
			Essential(:-, :y, 0),
			Essential(:-, :x, 0),
			Essential(:+, :x, 0),
			Essential(:+, :y, 0)
		),
		y = (
			Essential(:-, :y, 0),
			Essential(:-, :x, 0),
			Essential(:+, :x, 0),
			Essential(:+, :y, 0)
		)
	)
	bc_impl_mode = (
		x = (
			Essential(:-, :y,     f.x),
			Essential(:-, :x,     f.x),
			Essential(:+, :x,     f.x),
			Essential(:+, :y,     f.x)
		),
		y = (
			Essential(:-, :y,     f.y),
			Essential(:-, :x,     f.y),
			Essential(:+, :x,     f.y),
			Essential(:+, :y,     f.y)
		)
	)

	# helpers
	h   = Vector(p.n, motion_stags)
	r   = Vector(p.n, motion_stags)
	
	# modes
	m_1 = Vector(p.n, motion_stags)
	m_n = Vector(p.n, motion_stags)
	
	(λ_1, λ_n) = extremal_eigenmodes!(A, (m_1, m_n, bc_expl_mode), (r, bc_impl_mode), (h, f); bounds = bounds)
	
	plt1 = heatmap(x, y, m_1, "m_1", c = :davos)
	plt2 = heatmap(x, y, m_n, "m_n", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	readline()
	
	bc = (
		x = (
			Essential(:-, :y, -FD(f.x, :y)),
			Essential(:-, :x,     f.x),
			Essential(:+, :x,     f.x),
			Essential(:+, :y,     f.x-1)
		),
		y = (
			Essential(:-, :y,     f.y),
			Essential(:-, :x, -FD(f.y, :x)),
			Essential(:+, :x,  BD(f.y, :x)),
			Essential(:+, :y,     f.y)
		)
	)
	assign!(f, gen_rand(f), bounds)
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-5, maxit = 5000)

	display(plot(log10.(ε)))
	readline()

	assign!(b, (x = 1/(p.n[1]*p.n[2]), y = 0), bounds)
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-5, maxit = 5000)

	display(plot(log10.(ε)))
	readline()
	
	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
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
	
	   test_poisson(parameters(; parse_args(s)...))
	# test_elasticity(parameters(; parse_args(s)...))
	
	# "finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end