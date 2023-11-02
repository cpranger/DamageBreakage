#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

function test_poisson(p, ax_x, ax_y)
	# A(u) v = b
	u =  Field(p.n, div_stags)
	v = Vector(p.n, motion_stags)
	
	bc(u) = (
		"-y" => FD(  u, :y),
		"-x" =>   ( -u    ),
		"+x" =>   ( -u    ),
		"+y" =>   ( -u + 1)
	)

	f_u(u, v) = (divergence(v), bc(u))
	f_v(u, v) =  grad(u) - v

	sc = SchurComplement(f_u, f_v, v)

	# helpers
	w   = deepcopy(u)
	r   = deepcopy(u)
	h   = Tuple([deepcopy(u) for _ in 1:5])
		
	assign!(u, (0, bc(u)))
	
	newtonit!(sc, u, w, r, h; maxit = 30, atol = 1e-9)
	
	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, v, "v", c = :davos)
	plt3 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1, 3))
	display(plt)
end

function test_elastic(p, ax_x, ax_y)
	s(e) = (p.λ_0 * Tensor(MajorIdentity) + p.μ_0 * Tensor(MinorIdentity)) * e
	
	# A(u) v = b
	u = Vector(p.n, motion_stags)
	v = Tensor(p.n, Symmetric, strain_stags)
	
	bc(u) = (;
		x = (
			"-y" => FD(u.x, :y),
			"-x" =>   -u.x,
			"+x" =>   -u.x,
			"+y" =>  1-u.x
		),
		y = (
			"-y" =>   -u.y,
			"-x" =>   -u.y,
			"+x" =>   -u.y,
			"+y" =>   -u.y
		)
	)
	f_u(u, v) = (divergence(s(v)), bc(u))
	f_v(u, v) = symgrad(u) - v
	
	sc = SchurComplement(f_u, f_v, v)

	# helpers
	w   = deepcopy(u)
	r   = deepcopy(u)
	h   = Tuple([deepcopy(u) for _ in 1:5])
	
	assign!(u, (0, bc(u)))

	newtonit!(sc, u, w, r, h; maxit = 30, atol = 1e-9)

	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, v, "v", c = :davos)
	plt3 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1, 3))
	display(plt)
end


function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	h    =  1 / (n[1] - 2)
	
	μ_0  =  1
	λ_0  =  μ_0                      # assuming Poisson's ratio ν = 1/4
	k_0  =  λ_0 + (2/3)*μ_0          # bulk modulus

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
	
	assign!(x, fieldgen(i -> i/(p.n[1]-2) - 0.5))
	assign!(y, fieldgen(i -> i/(p.n[2]-2) - 0.5))
	
	# test_poisson(p, x, y)
	test_elastic(p, x, y)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end