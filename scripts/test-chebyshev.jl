#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

gen_ones() = fieldgen((_...) -> 1)
gen_rand() = fieldgen((_...) -> rand())

gen_ones(f::Field) = gen_ones()
gen_rand(f::Field) = gen_rand()

(gen_ones(v::Tensor{S, NamedTuple{N}}) where {S, N}) = (; zip(N, gen_ones() for n in N...)...)
(gen_rand(v::Tensor{S, NamedTuple{N}}) where {S, N}) = (; zip(N, gen_rand() for n in N...)...)

function extremal_eigenmodes!(A, (m_1, m_n, bc_expl), (r, bc_impl), (h, f))
	powerit_atol = 1e-3
	assign!(m_n, gen_rand(m_n))
	λ_n = powerit!(A, (m_n, m_1, bc_expl); maxit = 10000, atol = powerit_atol)
	λ_n *= 1 + 2 * powerit_atol # estimated to be larger  than actual λ_n
	λ_1 = λ_n / 10000 # estimated to be smaller than actual λ_1
	assign!(m_1, gen_ones(m_1))
	(λ_1, λ_n) = rayleighquotientit!(A, (m_1, m_n, bc_expl), (r, bc_impl), (h, f), (λ_1, λ_n); atol = 1e-6, maxit = 10000)
	return (λ_1, λ_n)
end

function test_poisson(p)
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i))
	assign!(y, fieldgen(i -> i))
	
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
	
	(λ_1, λ_n) = extremal_eigenmodes!(A, (m_1, m_n, bc_expl_mode), (r, bc_impl_mode), (h, f))
	
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
	
	assign!(f, gen_rand(f))
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, atol = 1e-5, maxit = 5000)

	display(plot(log10.(ε)))
	readline()

	assign!(b, 1/(p.n[1]*p.n[2]))
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, atol = 1e-5, maxit = 5000)

	display(plot(log10.(ε)))
	readline()

	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	readline()
end

function test_elasticity(p)
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i))
	assign!(y, fieldgen(i -> i))
	
	A    = x -> divergence(symgrad(x))

	# A f = b
	f   = Tensor(p.n, motion_stags)
	b   = Tensor(p.n, motion_stags)
	
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
	h   = Tensor(p.n, motion_stags)
	r   = Tensor(p.n, motion_stags)
	
	# modes
	m_1 = Tensor(p.n, motion_stags)
	m_n = Tensor(p.n, motion_stags)
	
	(λ_1, λ_n) = extremal_eigenmodes!(A, (m_1, m_n, bc_expl_mode), (r, bc_impl_mode), (h, f))
	
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
	assign!(f, gen_rand(f))
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, atol = 1e-5, maxit = 5000)

	display(plot(log10.(ε)))
	readline()

	assign!(b, (x = 1/(p.n[1]*p.n[2]), y = 0))
	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, atol = 1e-5, maxit = 5000)

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
	test_elasticity(parameters(; parse_args(s)...))
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end