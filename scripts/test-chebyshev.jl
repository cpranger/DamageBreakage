#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

function extremal_eigenmodes!(A, (m_1, m_n, bc_expl), (r, bc_impl), (h, f); bounds, gen_ones, gen_rand)
	powerit_atol = 1e-3
	assign!(m_n, gen_rand, bounds)
	λ_n = powerit!(A, (m_n, bc_expl), (m_1, bc_expl); bounds = bounds, maxit = *(bounds[2]...), atol = powerit_atol)
	λ_n *= 1 + 2 * powerit_atol # estimated to be larger  than actual λ_n
	λ_1 = λ_n / *(bounds[2]...) # estimated to be smaller than actual λ_1
	assign!(m_1, gen_ones, bounds)
	(λ_1, λ_n) = rayleighquotientit!(A, (m_1, m_n, bc_expl), (r, bc_impl), (h, f), (λ_1, λ_n); bounds = bounds, atol = 1e-6, maxit = max(bounds[2]...))
	return (λ_1, λ_n)
end

function test_poisson(x, y, p)
	bounds = (p.o, p.n)

	A = x -> divergence(grad(x))
	
	# A f = b
	f = Field(p.n, div_stags)
	b = Field(p.n, div_stags)

	bc_expl_mode = (
		Essential(:-, :y, 0),
		Essential(:-, :x, 0),
		Essential(:+, :x, 0),
		Essential(:+, :y, 0)
	)
	
	bc_impl_mode = (
		Essential(:-, :y, f),
		Essential(:-, :x, f),
		Essential(:+, :x, f),
		Essential(:+, :y, f)
	)

	# helpers
	h  = Field(p.n, div_stags)
	r  = Field(p.n, div_stags)
	
	# modes
	m_1 = Field(p.n, div_stags)
	m_n = Field(p.n, div_stags)
	
	gen_rand = fieldgen((_...) -> rand())
	gen_ones = 1

	(λ_1, λ_n) = extremal_eigenmodes!(A, (m_1, m_n, bc_expl_mode), (r, bc_impl_mode), (h, f); bounds = bounds, gen_ones = gen_ones, gen_rand = gen_rand)
	
	plt1 = heatmap(x, y, m_1, "m_1", c = :davos)
	plt2 = heatmap(x, y, m_n, "m_n", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	readline()

	bc = (
		Essential(:-, :y, f),
		Essential(:-, :x, f),
		Essential(:+, :x, f),
		Essential(:+, :y, f-1)
	)
	
	assign!(f, gen_rand, bounds)
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-5, maxit = 5000)

	display(plot(log10.(ε)))
	readline()

	assign!(b, 1/(p.n[1]*p.n[2]), bounds)
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-5, maxit = 5000)

	display(plot(log10.(ε)))
	readline()

	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	readline()
end

function test_elasticity(x, y, p)
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
	
	gen_rand = (
		x = fieldgen((_...) -> rand()),
		y = fieldgen((_...) -> rand())
	)
	gen_ones = (x = 1, y = 1)
	
	(λ_1, λ_n) = extremal_eigenmodes!(A, (m_1, m_n, bc_expl_mode), (r, bc_impl_mode), (h, f); bounds = bounds, gen_ones = gen_ones, gen_rand = gen_rand)
	
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
	assign!(f, gen_rand, bounds)
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
	
	p = parameters(; parse_args(s)...); pprintln(p)
	
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	test_poisson(x, y, p)
	test_elasticity(x, y, p)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end