#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane


function test__parabolic_scalar(p, ax)
	f(u) = divergence(grad(u)) / p.h^2
	bc(u) = (
		BC(-2, 0), BC(+2, 0), BC(-1, 0), BC(+1, 0)
	)
	bc_init(u) = (
		BC(-2, FD(u, :y)),
		BC(+2,  1-u     ),
		BC(-1,   -u     ),
		BC(+1,   -u     )
	) ./ p.h^2

	intg = tr_bdf2(x -> (f(x), bc(x)), Field(p.n, div_stags); h_t = .05)

	init!(intg, x -> (-x, bc_init(x)); newton_maxit = 30, newton_rtol = 1e-6)
	plt1 = heatmap(ax..., intg.y, "y", c = :davos)
	plt2 = heatmap(ax..., intg.e, "e", c = :davos)
	plt  = plot(plt1, plt2; layout = (1,2))
	display(plt)
	readline()

	for i in 1:20
		println("STEP $i")
		err = step!(intg; newton_maxit = 3, newton_rtol = 1e-6)
		println("time stepping error: $err.")
		plt1 = heatmap(ax..., intg.y, "u", c = :davos)
		plt2 = heatmap(ax..., intg.e, "e", c = :davos)
		plt  = plot(plt1, plt2; layout = (1,2))
		display(plt)
		readline()
	end
end

function test__parabolic_vector(p, ax)
	s(e)  = (p.λ_0 * Tensor(MajorIdentity) + p.μ_0 * Tensor(MinorIdentity)) * e
	f(u)  = divergence(s(symgrad(u))) / p.h^2
	bc(u) = (;
		x = (BC(-2, 0), BC(+2, 0), BC(-1, 0), BC(+1, 0)) ./ p.h^2,
		y = (BC(-2, 0), BC(+2, 0), BC(-1, 0), BC(+1, 0)) ./ p.h^2
	)
	bc_init(u) = (;
		x = (
			BC(-2,   -u.x),
			BC(+2,  1-u.x),
			BC(-1,   -u.x),
			BC(+1,   -u.x)
		) ./ p.h^2,
		y = (
			BC(-2,   -u.y),
			BC(+2,   -u.y),
			BC(-1,   -u.y),
			BC(+1,   -u.y)
		) ./ p.h^2
	)
	
	intg = tr_bdf2(x -> (f(x), bc(x)), Vector(p.n, motion_stags); h_t = .05)

	init!(intg, x -> (-x, bc_init(x)); newton_maxit = 30, newton_rtol = 1e-6)
	plt1 = heatmap(ax..., intg.y, "y", c = :davos)
	plt2 = heatmap(ax..., intg.e, "e", c = :davos)
	plt  = plot(plt1, plt2; layout = (1,2))
	display(plt)
	readline()

	for i in 1:20
		println("STEP $i")
		err = step!(intg; newton_maxit = 3, newton_rtol = 1e-6)
		println("time stepping error: $err.")
		plt1 = heatmap(ax..., intg.y, "u", c = :davos)
		plt2 = heatmap(ax..., intg.e, "e", c = :davos)
		plt  = plot(plt1, plt2; layout = (1,2))
		display(plt)
		readline()
	end
end

function test_hyperbolic_scalar(p, ax)
	u = Field(p.n, div_stags)
	v = Vector(p.n, motion_stags)
	
	bc      = u -> (BC(-2, FD(u, :y) / p.h), BC(+2, 0), BC(-1, 0), BC(+1, 0))
	bc_init = u -> (
		BC(-2, FD(u, :y)),
		BC(+2,  1-u     ),
		BC(-1,   -u     ),
		BC(+1,   -u     )
	) ./ p.h^2
	
	f((u, v),) = (
		(divergence(v) / p.h, bc(u)),
		       grad(u) / p.h
	)
	
	intg = tr_bdf2_schur(f, (u, v); h_t = .01)

	init!(intg, ((u, v),) -> ((-u, bc_init(u)), -v); newton_maxit = 30, newton_rtol = 1e-6)
	plt1 = heatmap(ax..., intg.y[1], "y", c = :davos)
	plt2 = heatmap(ax..., intg.e[1], "e", c = :davos)
	plt  = plot(plt1, plt2; layout = (1,2))
	display(plt)

	for i in 1:500
		println("STEP $i")
		err = step!(intg; newton_maxit = 3, newton_rtol = 1e-6)
		println("time stepping error: $err.")
		plt1 = heatmap(ax..., intg.y[1], "u", c = :davos)
		plt2 = heatmap(ax..., intg.e[1], "e", c = :davos)
		plt  = plot(plt1, plt2; layout = (1,2))
		display(plt)
	end
end

function test_hyperbolic_vector(p, ax)
	u = Vector(p.n, motion_stags)
	v = Tensor(p.n, Symmetric, strain_stags)
	
	w(v) = (p.λ_0 * Tensor(MajorIdentity) + p.μ_0 * Tensor(MinorIdentity)) * v
	
	bc = u -> (;
		x = (BC(-2, p.c_s * FD(u.x, :y) / p.h), BC(+2, 0), BC(-1, 0), BC(+1, 0)),
		y = (BC(-2, 0), BC(+2, 0), BC(-1, 0), BC(+1, 0))
	)
	bc_init = u -> (;
		x = (
			BC(-2, FD(u.x, :y)),
			BC(+2,  1-u.x     ),
			BC(-1,   -u.x     ),
			BC(+1,   -u.x     )
		) ./ p.h^2,
		y = (
			BC(-2,  -u.y),
			BC(+2,  -u.y),
			BC(-1,  -u.y),
			BC(+1,  -u.y)
		) ./ p.h^2
	)
	
	f((u, v),) = (
		(divergence(w(v)) / p.h, bc(u)),
		      symgrad(u)  / p.h
	)
	
	intg = tr_bdf2_schur(f, (u, v); h_t = .01)

	init!(intg, ((u, v),) -> ((-u, bc_init(u)), -v); newton_maxit = 30, newton_rtol = 1e-6)
	plt1 = heatmap(ax..., intg.y[1], "y", c = :davos)
	plt2 = heatmap(ax..., intg.e[1], "e", c = :davos)
	plt  = plot(plt1, plt2; layout = (1,2))
	display(plt)

	for i in 1:500
		println("STEP $i")
		err = step!(intg; newton_maxit = 3, newton_rtol = 1e-6)
		println("time stepping error: $err.")
		plt1 = heatmap(ax..., intg.y[1], "u", c = :davos)
		plt2 = heatmap(ax..., intg.e[1], "e", c = :davos)
		plt  = plot(plt1, plt2; layout = (1,2))
		display(plt)
	end
end

function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	h    =  1 / (n[1] - 2)
	l    =  (n .- 2) .* h            # physical domain size
	
	r_0  =  1
	μ_0  =  1
	λ_0  =  μ_0                      # assuming Poisson's ratio ν = 1/4
	k_0  =  λ_0 + (2/3)*μ_0          # Bulk modulus
	
	c_s  =  sqrt(μ_0 / r_0)          # solenoidal  wave speed (m/s)
	c_p  =  sqrt((λ_0+2*μ_0)/r_0)    # compressive wave speed (m/s)
	
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
	
	p = parameters(; parse_args(s)...); pprintln(p)
	
	# axes
	ax  = (
	    Field((p.n[1],), ((0,), (1,))),
	    Field((p.n[2],), ((0,), (1,)))
	)
	
	assign!(ax[1], fieldgen(i -> i*p.h - 0.5))
	assign!(ax[2], fieldgen(i -> i*p.h - 0.5))
	
	# test__parabolic_scalar(p, ax)
	# test__parabolic_vector(p, ax)
	# test_hyperbolic_scalar(p, ax)
	test_hyperbolic_vector(p, ax)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end