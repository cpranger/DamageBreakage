#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane
# using LinearAlgebra

function test_poisson(p, ax)
	u = Field(p.n, div_stags)
	h = Tuple([deepcopy(u) for _ in 1:5])
	b = deepcopy(u)
	
	f  = u -> divergence(grad(u))
	bc(u) = (
		BC(-2, FD(u, :y)),
		BC(+2,  1-u     ),
		BC(-1,   -u     ),
		BC(+1,   -u     ),
	)

	A = linearize(u -> (f(u), bc(u)), 0*u)
	assign!(b, -((f(0*u), bc(0*u))))
	
	assign!(u, (0, bc(0)))
	(λ, Λ, ε) = cg_pc_jacobi!(A, u, b; h = h, rtol = 1e-6, λtol = 1e-2, minit = 10, maxit = 100)
	
	r = deepcopy(u)
	assign!(r, b - A(u))

	plt1 = heatmap(ax..., u, "u", c = :davos)
	plt2 = heatmap(ax..., r, "r", c = :davos)
	display <| plot(plt1, plt2; layout = (1, 2))
	
	return
end

function test_elastic(p, ax)
	u   = Tensor(p.n, motion_stags)
	b   = deepcopy(u)
	h   = Tuple([deepcopy(u) for _ in 1:5])
	
	s(e) = (p.λ_0 * Tensor(MajorIdentity) + p.μ_0 * Tensor(MinorIdentity)) * e
	
	f  = u -> divergence(s(symgrad(u)))
	bc = u -> (;
		x = (
			BC(-2, FD(u.x, :y)),
			BC(+2,  1-u.x     ),
			BC(-1,   -u.x     ),
			BC(+1,   -u.x     )
		),
		y = (
			BC(-2,   -u.y     ),
			BC(+2,   -u.y     ),
			BC(-1,   -u.y     ),
			BC(+1,   -u.y     )
		)
	)
	
	A = linearize(u -> (f(u), bc(u)), 0*u)
	assign!(b, -((f(u), bc(u))))
	
	assign!(u, (0, bc(0*u)))
	(λ, Λ, ε) = cg_pc_jacobi!(A, u, b; h = h, rtol = 1e-6, λtol = 1e-2, minit = 10, maxit = 200)
	
	r = deepcopy(u)
	assign!(r, b - A(u))

	plt1 = heatmap(ax..., u, "u", c = :davos)
	plt2 = heatmap(ax..., r, "r", c = :davos)
	display <| plot(plt1, plt2; layout = (1, 2))
	
	return
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
	ax  = (
	    Field((p.n[1],), ((0,), (1,))),
	    Field((p.n[2],), ((0,), (1,)))
	)
	
	assign!(ax[1], fieldgen(i -> i*p.h - 0.5))
	assign!(ax[2], fieldgen(i -> i*p.h - 0.5))
	
	test_poisson(p, ax)
	# test_elastic(p, ax)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end