#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

# gen_ones() = fieldgen((_...) -> 1)
# gen_chss() = fieldgen((i, j) -> Int(mod(i,2) == mod(j,2)) - Int(mod(i,2) != mod(j,2)))
# gen_rand() = fieldgen((_...) -> rand())

# gen_ones(::Field) = gen_ones()
# gen_chss(::Field) = gen_chss()
# gen_rand(::Field) = gen_rand()

# (gen_ones(v::Tensor{S, NamedTuple{N, T}}) where {S, N, T}) = (; zip(N, [gen_ones() for n in N])...)
# (gen_chss(v::Tensor{S, NamedTuple{N, T}}) where {S, N, T}) = (; zip(N, [gen_chss() for n in N])...)
# (gen_rand(v::Tensor{S, NamedTuple{N, T}}) where {S, N, T}) = (; zip(N, [gen_rand() for n in N])...)

function test_poisson(p, ax)
	u     = Field(p.n, div_stags)
	f     = x -> divergence(grad(x))
	# f   = x -> divergence(abs(interpolate(x))*grad(x)) - p.h^2*x*(1 - x)
	bc(u) = (
		BC(-2, FD(u, :y)),
		BC(+2,  1-u     ),
		BC(-1,   -u     ),
		BC(+1,   -u     ),
	)
	
	# helpers
	r    = deepcopy(u)
	h    = Tuple([deepcopy(u) for _ in 1:6])
		
	assign!(u, (0, bc(0)))
	newtonit!(x -> (f(x), bc(x)), u, r, h; maxit = 30, rtol = 1e-9)
	
	plt1 = heatmap(ax..., u, "u", c = :davos)
	plt2 = heatmap(ax..., r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	display(plt)
end

function test_elastic(p, ax)
	u    = Tensor(p.n, motion_stags)
	s(e) = (p.λ_0 * Tensor(MajorIdentity) + p.μ_0 * Tensor(MinorIdentity)) * e
	f    = u -> divergence(s(symgrad(u)))
	bc   = u -> (;
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
	
	# helpers
	r   = deepcopy(u)
	h   = Tuple([deepcopy(u) for _ in 1:6])
		
	assign!(u, (0, bc(0*u)))
	newtonit!(u -> (f(u), bc(u)), u, r, h; maxit = 30, rtol = 1e-9)

	plt1 = heatmap(ax..., u, "u", c = :davos)
	plt2 = heatmap(ax..., r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
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