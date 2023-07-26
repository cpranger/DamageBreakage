#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

gen_ones() = fieldgen((_...) -> 1)
gen_chss() = fieldgen((i, j) -> Int(mod(i,2) == mod(j,2)) - Int(mod(i,2) != mod(j,2)))
gen_rand() = fieldgen((_...) -> rand())

gen_ones(::Field) = gen_ones()
gen_chss(::Field) = gen_chss()
gen_rand(::Field) = gen_rand()

(gen_ones(v::Tensor{S, NamedTuple{N, T}}) where {S, N, T}) = (; zip(N, [gen_ones() for n in N])...)
(gen_chss(v::Tensor{S, NamedTuple{N, T}}) where {S, N, T}) = (; zip(N, [gen_chss() for n in N])...)
(gen_rand(v::Tensor{S, NamedTuple{N, T}}) where {S, N, T}) = (; zip(N, [gen_rand() for n in N])...)

linearize(f, x; h = eps(Float32)) = v -> imag(f(x + h * im * v)) / h

function test_poisson(p)
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	bounds = (p.o, p.n)

	# A(u) v = b
	u = Field(p.n, div_stags)
	v = Field(p.n, div_stags)
	b = Field(p.n, div_stags)

	F(x) = divergence(abs(interpolate(x))*grad(x)) - p.h^2*x*(1 - x)
	
	bc_expl_mode = (
		Essential(:-, :y, 0),
		Essential(:-, :x, 0),
		Essential(:+, :x, 0),
		Essential(:+, :y, 0)
	)
	
	bc_impl_mode = (
		Essential(:-, :y, v),
		Essential(:-, :x, v),
		Essential(:+, :x, v),
		Essential(:+, :y, v)
	)

	# helpers
	h  = Field(p.n, div_stags)
	f  = Field(p.n, div_stags)
	r  = Field(p.n, div_stags)
	
	# modes
	m_1 = Field(p.n, div_stags)
	m_n = Field(p.n, div_stags)
	
	bc_expl = (
		Essential(:-, :y, 0),
		Essential(:-, :x, 0),
		Essential(:+, :x, 0),
		Essential(:+, :y, 1)
	)
	
	bc_impl = (
		Essential(:-, :y, v),
		Essential(:-, :x, v),
		Essential(:+, :x, v),
		Essential(:+, :y, v)
	)
	
	assign!((u, bc_expl), gen_ones(u), bounds)
	
	assign!(m_n, gen_chss(m_n), bounds)
	assign!(m_1, gen_ones(m_1), bounds)
	
	λ_n  = powerit!(linearize(F, u), 0,   (m_n, h, bc_expl_mode); bounds = bounds, maxit = *(bounds[2]...), atol = 1e-3)
	λ_1  = powerit!(linearize(F, u), λ_n, (m_1, h, bc_expl_mode); bounds = bounds, maxit = *(bounds[2]...), atol = 1e-3)
	
	Meta.@show (λ_1, λ_n)
	
	λ_n = rayleighquotientit!(linearize(F, u), λ_n + λ_1, (m_n, bc_expl_mode), (r, bc_impl_mode), (h, f), (-λ_1, -λ_n); bounds = bounds, atol = 1e-6, maxit = 10, chebymaxit = max(bounds[2]...))
	λ_1 = rayleighquotientit!(linearize(F, u), 0,         (m_1, bc_expl_mode), (r, bc_impl_mode), (h, f), (+λ_1, +λ_n); bounds = bounds, atol = 1e-6, maxit = 10, chebymaxit = max(bounds[2]...))
	
	Meta.@show (λ_1, λ_n)
	
	plt1 = heatmap(x, y, m_1, "m_1", c = :davos)
	plt2 = heatmap(x, y, m_n, "m_n", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	
	assign!(v, m_1, bounds)
	newtonit!(F, u, (v, bc_expl_mode), (r, bc_impl), ((λ_1, m_1), (λ_n, m_n)), (h, f); bounds = bounds, atol = 1e-6, maxit = 20)

	# plt1 = heatmap(x, y, u, "u", c = :davos)
	# plt2 = heatmap(x, y, r, "r", c = :davos)
	# plt  = plot(plt1, plt2; layout = (1, 2))
	# display(plt)
end

# function test_elasticity(p)
# 	# axes
# 	x  = Field((p.n[1],), ((0,), (1,)))
# 	y  = Field((p.n[2],), ((0,), (1,)))
	
# 	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
# 	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
# 	bounds = (p.o, p.n)

# 	A    = x -> divergence(symgrad(x))

# 	# A f = b
# 	f   = Vector(p.n, motion_stags)
# 	b   = Vector(p.n, motion_stags)
	
# 	bc_expl_mode = (
# 		x = (
# 			Essential(:-, :y, 0),
# 			Essential(:-, :x, 0),
# 			Essential(:+, :x, 0),
# 			Essential(:+, :y, 0)
# 		),
# 		y = (
# 			Essential(:-, :y, 0),
# 			Essential(:-, :x, 0),
# 			Essential(:+, :x, 0),
# 			Essential(:+, :y, 0)
# 		)
# 	)
# 	bc_impl_mode = (
# 		x = (
# 			Essential(:-, :y,     f.x),
# 			Essential(:-, :x,     f.x),
# 			Essential(:+, :x,     f.x),
# 			Essential(:+, :y,     f.x)
# 		),
# 		y = (
# 			Essential(:-, :y,     f.y),
# 			Essential(:-, :x,     f.y),
# 			Essential(:+, :x,     f.y),
# 			Essential(:+, :y,     f.y)
# 		)
# 	)

# 	# helpers
# 	h   = Vector(p.n, motion_stags)
# 	r   = Vector(p.n, motion_stags)
	
# 	# modes
# 	m_1 = Vector(p.n, motion_stags)
# 	m_n = Vector(p.n, motion_stags)
	
# 	(λ_1, λ_n) = extremal_eigenmodes!(A, (m_1, m_n, bc_expl_mode), (r, bc_impl_mode), (h, f); bounds = bounds, atol = 1e-5)
	
# 	plt1 = heatmap(x, y, m_1, "m_1", c = :davos)
# 	plt2 = heatmap(x, y, m_n, "m_n", c = :davos)
# 	plt  = plot(plt1, plt2; layout = (1, 2))
# 	display(plt)
# 	readline()
	
# 	bc = (
# 		x = (
# 			Essential(:-, :y, -FD(f.x, :y)),
# 			Essential(:-, :x,     f.x),
# 			Essential(:+, :x,     f.x),
# 			Essential(:+, :y,     f.x-1)
# 		),
# 		y = (
# 			Essential(:-, :y,     f.y),
# 			Essential(:-, :x, -FD(f.y, :x)),
# 			Essential(:+, :x,  BD(f.y, :x)),
# 			Essential(:+, :y,     f.y)
# 		)
# 	)
# 	assign!(f, gen_rand(f), bounds)
# 	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-5, maxit = 5000)

# 	display(plot(log10.(ε)))
# 	readline()

# 	assign!(b, (x = 1/(p.n[1]*p.n[2]), y = 0), bounds)
# 	ε = chebyshev!(A, f, b, (r, bc); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-5, maxit = 5000)

# 	display(plot(log10.(ε)))
# 	readline()
	
# 	plt1 = heatmap(x, y, f, "f", c = :davos)
# 	plt2 = heatmap(x, y, r, "r", c = :davos)
# 	plt  = plot(plt1, plt2; layout = (1, 2))
	
# 	display(plt)
	
# 	readline()
# end

function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	o    =  n .- n .+ 1              # logical origin
	h    =  1 / (n[1] - 2)
	
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
	#    readline()
	# test_elasticity(parameters(; parse_args(s)...))
	# readline()
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end