#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane


function test__parabolic_scalar(p, ax_x, ax_y)
	
	u =  Field(p.n, div_stags)
	
	f(u) = -divergence(grad(u))

	bc(u0 = 0) = u -> (
		"-y" => -FD(u,  :y),
		"-x" =>     u,
		"+x" =>     u,
		"+y" =>     u - u0
	)

	intg = tr_bdf2(x -> (f(x), bc(0)(x)), u; h_t = 1/p.h^2)

	assign!(u, (0, -bc(1)(0)))

	newtonit!(x -> (-f(x), -bc(1)(x)), u, intg.dw, intg.r, intg.h; maxit = 30, atol = 1e-5)
	
	display(heatmap(ax_x, ax_y, u, "u", c = :davos))

	for i in 1:1
		println("STEP $i")
		err = step!(intg; newton_maxit = 3, newton_atol = 1e-9)
		println("time stepping error: $err.")
		plt1 = heatmap(ax_x, ax_y,      u, "u", c = :davos)
		plt2 = heatmap(ax_x, ax_y, intg.e, "e", c = :davos)
		# plt3 = heatmap(ax_x, ax_y, intg.r, "r", c = :davos)
		plt  = plot(plt1, plt2#=, plt3=#; layout = (1,2))
		display(plt)
	end
	
end

function test__parabolic_vector(p, ax_x, ax_y)
	
	v = Vector(p.n, motion_stags)
	
	s(e) = (p.λ_0 * MajorIdentity() + p.μ_0 * MinorIdentity()) * e

	f(v) = divergence(s(symgrad(v))) / p.h^2

	bc(v) = (;
		x = (
			"-y" => FD(  v.x, :y) / p.h^2,
			"-x" =>   ( -v.x    ) / p.h^2,
			"+x" =>   ( -v.x    ) / p.h^2,
			"+y" =>   (1-v.x    ) / p.h^2
		),
		y = (
			"-y" =>   ( -v.y    ) / p.h^2,
			"-x" =>   ( -v.y    ) / p.h^2,
			"+x" =>   ( -v.y    ) / p.h^2,
			"+y" =>   ( -v.y    ) / p.h^2
		)
	)

	intg = tr_bdf2(f, u)
	
	err = step!(intg; bc = bc, newton_maxit = 30, newton_atol = 1e-9)

	println("time stepping error: $err.")

	plt1 = heatmap(ax_x, ax_y,      v, "v", c = :davos)
	plt2 = heatmap(ax_x, ax_y, intg.e, "e", c = :davos)
	plt3 = heatmap(ax_x, ax_y, intg.r, "r", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1, 3))
	display(plt)
end

function test_hyperbolic_scalar(p, ax_x, ax_y)
	
	u =  Field(p.n, div_stags)
	v = Vector(p.n, motion_stags)
	
	f_u(u, v) = p(divergence(v)) / p.h
	f_v(u, v) = grad(u) / p.h - v
	bc_(u)    = (
		"-y" => FD( u, :y) / p.h^2,
		"-x" =>   (-u    ) / p.h^2,
		"+x" =>   (-u    ) / p.h^2,
		"+y" =>   (-u + 1) / p.h^2
	)

	s = SchurComplement(f_u, f_v, v)

	intg = tr_bdf2(s, u)
	
	err = step!(intg; bc = bc_, newton_maxit = 30, newton_atol = 1e-9)

	println("time stepping error: $err.")

	plt1 = heatmap(ax_x, ax_y,      u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, intg.e, "e", c = :davos)
	plt3 = heatmap(ax_x, ax_y, intg.r, "r", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1, 3))
	display(plt)
end

function test_hyperbolic_vector(p, ax_x, ax_y)
	
	v = Vector(p.n, motion_stags)
	e = Tensor(p.n, Symmetric, strain_stags)
	
	s(e) = (p.λ_0 * MajorIdentity() + p.μ_0 * MinorIdentity()) * e

	f_v(v, e) = divergence(s(e)) / p.h
	f_e(v, e) = symgrad(v) / p.h - e
	bc_(v)    = (;
		x = (
			"-y" => FD(  v.x, :y) / p.h^2,
			"-x" =>   ( -v.x    ) / p.h^2,
			"+x" =>   ( -v.x    ) / p.h^2,
			"+y" =>   (1-v.x    ) / p.h^2
		),
		y = (
			"-y" =>   ( -v.y    ) / p.h^2,
			"-x" =>   ( -v.y    ) / p.h^2,
			"+x" =>   ( -v.y    ) / p.h^2,
			"+y" =>   ( -v.y    ) / p.h^2
		)
	)
	
	s = SchurComplement(f_v, f_e, e)

	intg = tr_bdf2(s, u)
	
	err = step!(intg; bc = bc_, newton_maxit = 30, newton_atol = 1e-9)

	println("time stepping error: $err.")

	plt1 = heatmap(ax_x, ax_y, v, "v", c = :davos)
	plt2 = heatmap(ax_x, ax_y, e, "e", c = :davos)
	plt3 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1, 3))
	display(plt)
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
	ax_x  = Field((p.n[1],), ((0,), (1,)))
	ax_y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(ax_x, fieldgen(i -> i))
	assign!(ax_y, fieldgen(i -> i))
	
	test__parabolic_scalar(p, ax_x, ax_y)
	# test__parabolic_vector(p, ax_x, ax_y)
	# test_hyperbolic_scalar(p, ax_x, ax_y)
	# test_hyperbolic_vector(p, ax_x, ax_y)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end