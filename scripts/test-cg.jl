#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

function test_poisson(p, ax_x, ax_y)
	u   = Field(p.n, div_stags)
	h_1 = deepcopy(u)
	h_2 = deepcopy(u)
	r   = deepcopy(u)
	b   = deepcopy(u)
	pc  = deepcopy(u)
	
	f  = u -> divergence(grad(u)) - 1/1000
	bc = u -> (
		"-y" => FD(  u, :y),
		"-x" =>   ( -u    ),
		"+x" =>   ( -u    ),
		"+y" =>   ( -u+1  )
	)

	A = linearize(u -> (f(u), bc(u)), 0*u)
	assign!(b, (f(0*u), bc(0*u)))
	assign!(b, -b)

	assign!(u, (0, bc(0*u)))
	(λ, Λ, ε) = cg!(A, u, b; r = r, p = h_1, q = h_2, rtol = 1e-8, λtol = 1e-2, minit = 100, maxit = 1000)
	# (λ, Λ, ε) = cg_pc_jacobi!(A, u, b; d = pc, r = r, p = h_1, q = h_2, rtol = 1e-8, λtol = 1e-2, minit = 100, maxit = 1000)

	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	display <| plot(plt1, plt2; layout = (1, 2))
	
	readline()
	return
end

function test_elastic(p, ax_x, ax_y)
	u   = Vector(p.n, motion_stags)
	h_1 = deepcopy(u)
	h_2 = deepcopy(u)
	r   = deepcopy(u)
	b   = deepcopy(u)
	
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
	
	A = linearize(u -> (f(u), bc(u)), 0*u)
	assign!(b, (f(0*u), bc(0*u)))
	assign!(b, -b)
	
	assign!(u, (0, bc(0*u)))
	(λ, Λ, ε) = cg!(A, u, b; r = r, p = h_1, q = h_2, rtol = 1e-8, λtol = 1e-2, minit = 100, maxit = 1000)
	
	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	display <| plot(plt1, plt2; layout = (1, 2))
	
	readline()
	return
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
	
	p = parameters(; parse_args(s)...)

	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i))
	assign!(y, fieldgen(i -> i))
	
	test_poisson(p, x, y)
	# test_elastic(p, x, y)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end