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

linearize(F, x; h = eps(Float32)) = v -> imag(F(x + h * im * v)) / h

linearize(F, a::Tuple, x, b::Tuple; h = eps(Float32)) = v -> imag(F(a..., x + h * im * v, b...)) / h

function test_poisson(p)
	# # axes
	ax_x  = Field((p.n[1],), ((0,), (1,)))
	ax_y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(ax_x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(ax_y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	bounds = (p.o, p.n)

	x  = Field(p.n, div_stags)
	dx = Field(p.n, div_stags)
	u  = Field(p.n, div_stags)

	y  = Vector(p.n, motion_stags)
	dy = Vector(p.n, motion_stags)
	v  = Vector(p.n, motion_stags)
	
	# |A B| |x| _ |u|
	# |C D| |y| ‾ |v|

	# S(x) := A(x) - B(1/D(C(x)))
	# S(x)  = u - 1/D(C(v))
	# y     = 1/D(v - C(x))
	
	Q(x) = divergence(abs(interpolate(x))*grad(x))# - p.h^2*x*(1 - x)
	
	AB(x, y) = #=p.h^2*x*(x - 1) +=# divergence(y)
	CD(x, y) = abs(interpolate(x))*grad(x) - y

	A = linearize(AB, (  ), x, (y,)) # == diag
	B = linearize(AB, (x,), y, (  ))
	C = linearize(CD, (  ), x, (y,))
	D = linearize(CD, (x,), y, (  )) # == diag

	S(x)     = A(x) - B(1/D(1/C(x)))
	Y(x, v)  = 1/D(v - C(x))
	
	bc_expl_mode = (
		Essential(:-, :y, 0),
		Essential(:-, :x, 0),
		Essential(:+, :x, 0),
		Essential(:+, :y, 0)
	)
	
	bc_impl_mode = (
		Essential(:-, :y, dx),
		Essential(:-, :x, dx),
		Essential(:+, :x, dx),
		Essential(:+, :y, dx)
	)

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

	assign!(mx_1,  Mode(mx_1, bc_x)[+1, +1].gen, bounds)
	# assign!(mx_n,  Mode(mx_n, bc_x)[-1, -1].gen, bounds)
	assign!(my_1, pMode(my_1, bc_y)[+1, +1].gen, bounds)
	# assign!(my_n, pMode(my_n, bc_y)[-1, -1].gen, bounds)
	
	linearized = true

	if linearized
		A = linearize(AB, (  ), mx_1, (my_1,))
		assign!(( x, bc_x), A(mx_1), bounds)
	else
		assign!(( x, bc_x), 0*mx_1, bounds)
	end

	if linearized
		B = linearize(AB, (mx_1,), my_1, (  ))
		assign!((dx, bc_x), B(my_1), bounds)
	else
		assign!((dx, bc_x), divergence(my_1), bounds)
	end

	if linearized
		C = linearize(CD, (  ), mx_1, (my_1,))
		assign!(( y, bc_y), C(mx_1), bounds)
	else
		assign!(( y, bc_y), abs(interpolate(mx_1))*grad(mx_1), bounds)
	end

	if linearized
		D = linearize(CD, (mx_1,), my_1, (  ))
		assign!((dy, bc_y), D(my_1), bounds)
	else
		assign!((dy, bc_y),  -my_1,  bounds)
	end
	
	plt1 = heatmap(ax_x, ax_y,   x,    "x"  , c = :davos)
	plt2 = heatmap(ax_x, ax_y,  dx,   "dx"  , c = :davos)
	plt3 = heatmap(ax_x, ax_y,   y.x,  "y.x", c = :davos)
	plt4 = heatmap(ax_x, ax_y,   y.y,  "y.y", c = :davos)
	plt5 = heatmap(ax_x, ax_y,  dy.x, "dy.x", c = :davos)
	plt6 = heatmap(ax_x, ax_y,  dy.y, "dy.y", c = :davos)
	plt  = plot(plt1, plt2, plt3, plt4, plt5, plt6; layout = (3, 2))
	display(plt)
	readline()

	if linearized
		A = linearize(AB, (  ), mx_1, (my_1,)) # == diag
		B = linearize(AB, (mx_1,), my_1, (  ))
		C = linearize(CD, (  ), mx_1, (my_1,))
		D = linearize(CD, (mx_1,), my_1, (  )) # == diag

		# S(x)     = A(x) - B(1/D(1/C(x)))
		assign!(( x, Natural()), S(mx_1), bounds)
	else
		assign!(( x, Natural()), Q(mx_1), bounds)
	end
	
	plt1 = heatmap(ax_x, ax_y,   x,    "x"  , c = :davos)
	display(plt1)

	λ_n  = powerit!(S, 0,   (mx_n, bc_expl_mode), (h, bc_expl_mode); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
	λ_1  = powerit!(S, λ_n, (mx_1, bc_expl_mode), (h, bc_expl_mode); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
	
	Meta.@show (λ_1, λ_n)
	
	plt1 = heatmap(ax_x, ax_y, mx_1, "m_1", c = :davos)
	plt2 = heatmap(ax_x, ax_y, mx_n, "m_n", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
	return; readline()
	
	λ_n = rayleighquotientit!(S, λ_n + λ_1, (mx_n, bc_expl_mode), (r, bc_impl_mode), (h, f), (-λ_1, -λ_n); bounds = bounds, atol = 1e-6, maxit = 10, chebymaxit = max(bounds[2]...))
	λ_1 = rayleighquotientit!(S, 0,         (mx_1, bc_expl_mode), (r, bc_impl_mode), (h, f), (+λ_1, +λ_n); bounds = bounds, atol = 1e-6, maxit = 10, chebymaxit = max(bounds[2]...))
	
	Meta.@show (λ_1, λ_n)
	
	# plt1 = heatmap(x, y, m_1, "m_1", c = :davos)
	# plt2 = heatmap(x, y, m_n, "m_n", c = :davos)
	# plt  = plot(plt1, plt2; layout = (1, 2))
	# display(plt)
	
	# assign!(v, m_1, bounds)
	chebyshev!(S, dx, u, (r, bc_impl); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = 1e-6, maxit = 10*max(bounds[2]...))
	# newtonit!(F, u, b, (v, bc_expl_mode), (r, bc_impl), ((λ_1, m_1), (λ_n, m_n)), (h, f); bounds = bounds, atol = 1e-6, maxit = 20)

	# plt1 = heatmap(x, y, u, "u", c = :davos)
	# plt2 = heatmap(x, y, r, "r", c = :davos)
	# plt  = plot(plt1, plt2; layout = (1, 2))
	# display(plt)
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