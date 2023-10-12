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

function test_poisson(p, ax_x, ax_y)
	# A(u) v = b
	u = Field(p.n, div_stags)
	v = Field(p.n, div_stags)

	f  = x -> divergence(grad(x))
	# f  = x -> divergence(abs(interpolate(x))*grad(x)) - p.h^2*x*(1 - x)
	bc = x -> (
		"-y" => FD(  x, :y),
		"-x" =>   ( -x    ),
		"+x" =>   ( -x    ),
		"+y" =>   ( -x+1  )
	)
	
	# helpers
	r   = Field(p.n, div_stags)
	h_1 = Field(p.n, div_stags)
	h_2 = Field(p.n, div_stags)
	h_3 = Field(p.n, div_stags)
		
	assign!(u, 1)
	
	newtonit!(x -> (f(x), bc(x)), u, v, r, (h_1, h_2, h_3); maxit = 30, atol = 1e-9)

	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
end

function test_elastic(p, ax_x, ax_y)
	# A(u) v = b
	u = Vector(p.n, motion_stags)
	v = Vector(p.n, motion_stags)

	f  = u -> divergence(symgrad(u))
	bc = u -> (;
		x = (
			"-y" => FD(  u.x, :y),
			"-x" =>   ( -u.x    ),
			"+x" =>   ( -u.x    ),
			"+y" =>   (1-u.x    )
		),
		y = (
			"-y" =>   ( -u.y    ),
			"-x" =>   ( -u.y    ),
			"+x" =>   ( -u.y    ),
			"+y" =>   ( -u.y    )
		)
	)
	
	# helpers
	r   = Vector(p.n, motion_stags)
	h_1 = Vector(p.n, motion_stags)
	h_2 = Vector(p.n, motion_stags)
	h_3 = Vector(p.n, motion_stags)
		
	assign!(u, gen_ones(u))
	
	newtonit!(u -> (f(u), bc(u)), u, v, r, (h_1, h_2, h_3); maxit = 30, atol = 1e-9)

	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
end

function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
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
	
	p = parameters(; parse_args(s)...)

	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i))
	assign!(y, fieldgen(i -> i))
	
	test_poisson(p, x, y)
	test_elastic(p, x, y)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end