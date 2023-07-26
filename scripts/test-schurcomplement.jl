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

function test_poisson(p)
	# # axes
	ax_x  = Field((p.n[1],), ((0,), (1,)))
	ax_y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(ax_x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(ax_y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	bounds = (p.o, p.n)

	x    = Field(p.n, div_stags)
	dx   = Field(p.n, div_stags)
	x_h  = Field(p.n, div_stags)
	
	y    = Vector(p.n, motion_stags)
	y_r  = Vector(p.n, motion_stags)
	y_h  = Vector(p.n, motion_stags)
	
	q(x) = abs(interpolate(x))*grad(x)
	
	# F(x) = divergence(q(x))# - p.h^2*x*(1 - x)
	
	F_x(x, y) = #=p.h^2*x*(x - 1) +=# divergence(y)
	F_y(x, y) = q(x) - y

	S = SchurComplement(F_x, F_y, y, y_r)
	
	# helpers
	h  = Field(p.n, div_stags)
	f  = Field(p.n, div_stags)
	r  = Field(p.n, div_stags)
	
	# modes
	mx_1 = Field(p.n, div_stags)
	mx_n = Field(p.n, div_stags)

	my_1 = Vector(p.n, motion_stags)
	my_n = Vector(p.n, motion_stags)
	
	bc_expl = (
		Essential(:-, :y, 0),
		Essential(:-, :x, 0),
		Essential(:+, :x, 0),
		Essential(:+, :y, 1)
	)
	
	bc_impl = (
		Essential(:-, :y, dx),
		Essential(:-, :x, dx),
		Essential(:+, :x, dx),
		Essential(:+, :y, dx)
	)
	
	assign!((x, bc_expl), gen_ones(x), bounds)
	
	bc_x = Essential()
	bc_y = ImpermeableFreeSlip().v
	
	assign!(mx_1, gen_rand(mx_1)#=Mode(mx_1, bc_x)[+1, +1].gen=#, bounds)
	assign!(mx_n, gen_rand(mx_n)#=Mode(mx_n, bc_x)[-1, -1].gen=#, bounds)
	# assign!(my_1, pMode(my_1, bc_y)[+1, +1].gen, bounds)
	# assign!(my_n, pMode(my_n, bc_y)[-1, -1].gen, bounds)
	
	# plt1 = heatmap(ax_x, ax_y, mx_1, "mx_1", c = :davos)
	# plt2 = heatmap(ax_x, ax_y, mx_n, "mx_n", c = :davos)
	# plt = plot(plt1, plt2, layout = (1,2))
	# display(plt)

	bc_expl_mode = (
		Essential(:-, :y, 0),
		Essential(:-, :x, 0),
		Essential(:+, :x, 0),
		Essential(:+, :y, 0),
	)
	
	bc_impl_mode = (
		Essential(:-, :y, dx),
		Essential(:-, :x, dx),
		Essential(:+, :x, dx),
		Essential(:+, :y, dx),
	)

	λ_n  = powerit!(linearize(S, x), 0,   (mx_n, x_h, bc_expl_mode); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
	λ_1  = powerit!(linearize(S, x), λ_n, (mx_1, x_h, bc_expl_mode); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
	
	Meta.@show (λ_1, λ_n)

	plt1 = heatmap(ax_x, ax_y, mx_1, "m_1", c = :davos)
	plt2 = heatmap(ax_x, ax_y, mx_n, "m_n", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	
	λ_n = rayleighquotientit!(linearize(S, x), λ_n + λ_1, (mx_n, bc_expl_mode), (r, bc_impl_mode), (h, f), (-λ_1, -λ_n); bounds = bounds, atol = 1e-6, maxit = 10, chebymaxit = max(bounds[2]...))
	λ_1 = rayleighquotientit!(linearize(S, x), 0,         (mx_1, bc_expl_mode), (r, bc_impl_mode), (h, f), (+λ_1, +λ_n); bounds = bounds, atol = 1e-6, maxit = 10, chebymaxit = max(bounds[2]...))
	
	Meta.@show (λ_1, λ_n)
	
	plt1 = heatmap(ax_x, ax_y, mx_1, "m_1", c = :davos)
	plt2 = heatmap(ax_x, ax_y, mx_n, "m_n", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	
	# assign!(v, m_1, bounds)
	# chebyshev!(S, dx, u, (r, bc_impl); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-6, maxit = 10*max(bounds[2]...))
	newtonit!(S, x, (dx, bc_expl_mode), (r, bc_impl), ((λ_1, mx_1), (λ_n, mx_n)), (h, f); bounds = bounds, atol = 1e-6, maxit = 7)

	plt1 = heatmap(ax_x, ax_y, x, "x", c = :davos)
	plt2 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
end


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