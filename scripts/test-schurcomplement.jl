#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

function test_poisson(p)
	# axes
	ax_x  = Field((p.n[1],), ((0,), (1,)))
	ax_y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(ax_x, fieldgen(i -> i))
	assign!(ax_y, fieldgen(i -> i))
	
	# A(u) v = b
	u =  Field(p.n, div_stags)
	v = Vector(p.n, motion_stags)
	
	f_u(u, v) = (
		divergence(v),
		(
			Essential(:-, :y, FD(u, :y)),
			Essential(:-, :x, -u),
			Essential(:+, :x, -u),
			Essential(:+, :y, -u + 1)
		)
	)

	f_v(u, v) = grad(u) - v

	s = SchurComplement(f_u, f_v, v)

	# helpers
	w   = Field(p.n, div_stags)
	r   = Field(p.n, div_stags)
	h_1 = Field(p.n, div_stags)
	h_2 = Field(p.n, div_stags)
	h_3 = Field(p.n, div_stags)
		
	assign!(u, 1)
	
	newtonit!(s, u, w, r, (h_1, h_2, h_3); maxit = 30, atol = 1e-9)

	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, v, "v", c = :davos)
	plt3 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1, 3))
	display(plt)
end

function test_elastic(p)
	# # axes
	ax_x  = Field((p.n[1],), ((0,), (1,)))
	ax_y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(ax_x, fieldgen(i -> i))
	assign!(ax_y, fieldgen(i -> i))
	
	# A(u) v = b
	u = Vector(p.n, motion_stags)
	v = Tensor(p.n, Symmetric, strain_stags)
	
	f_u(u, v) = (
		divergence(v),
		(;
			x = (
				Essential(:-, :y, FD(u.x, :y)),
				Essential(:-, :x,   -u.x),
				Essential(:+, :x,   -u.x),
				Essential(:+, :y,  1-u.x)
			),
			y = (
				Essential(:-, :y,   -u.y),
				Essential(:-, :x,   -u.y),
				Essential(:+, :x,   -u.y),
				Essential(:+, :y,   -u.y)
			)
		)
	)

	f_v(u, v) = symgrad(u) - v
	
	s = SchurComplement(f_u, f_v, v)

	# helpers
	w   = Vector(p.n, motion_stags)
	r   = Vector(p.n, motion_stags)
	h_1 = Vector(p.n, motion_stags)
	h_2 = Vector(p.n, motion_stags)
	h_3 = Vector(p.n, motion_stags)
		
	assign!(u, 1)
	
	newtonit!(s, u, w, r, (h_1, h_2, h_3); maxit = 30, atol = 1e-9)

	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, v, "v", c = :davos)
	plt3 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1, 3))
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
	
	# test_poisson(parameters(; parse_args(s)...))
	#    readline()
	test_elastic(parameters(; parse_args(s)...))
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end