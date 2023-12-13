#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

function display_2d(ax, intg::tr_bdf2)
	plt1 = heatmap(ax..., intg.y, "y", c = :davos)
	plt2 = heatmap(ax..., intg.e, "e", c = :davos)
	plt  = plot(plt1, plt2; layout = (1,2))
	display(plt)
end

function display_2d(ax, intg::tr_bdf2_schur)
	plt1 = heatmap(ax..., intg.y[1], "y", c = :davos)
	plt2 = heatmap(ax..., intg.e[1], "e", c = :davos)
	plt  = plot(plt1, plt2; layout = (1,2))
	display(plt)
end

lid_bc(s #=true for static=#, c #=absorption coeff=#) = u -> (
	BC(-2,  s*( -u)#=c*FD(u, :y)=#),
	BC(+2,  s*(s-u)),
	BC(-1,  s*( -u)),
	BC(+1,  s*( -u))
)

zero_bc(s #=true for static=#) = u -> (
	BC(-2, s*( -u)),
	BC(+2, s*( -u)),
	BC(-1, s*( -u)),
	BC(+1, s*( -u))
)

function lid_driven(rank, bolicity, imex, p, ax; duration = 4, rtol = 1e-4, atol = 1e-4)
	if     rank == :scalar
		u =  Field(p.n, div_stags)
		v = Vector(p.n, motion_stags)
		
		bc = (u; static = false) -> lid_bc(static, 1/p.h)(u)
		
		w = v -> v
		_grad = grad
	
	
	elseif rank == :vector
		u = Vector(p.n, motion_stags)
		v = Tensor(p.n, Symmetric, strain_stags)
		
		bc = (u; static = false) -> (;
			x =  lid_bc(static, p.c_s/p.h)(u.x),
			y = zero_bc(static           )(u.y)
		)

		w = v -> (p.λ_0 * Tensor(MajorIdentity) + p.μ_0 * Tensor(MinorIdentity)) * v
		_grad = symgrad
	else
		error("argument 'rank' is neither of (:scalar, :vector)")
	end
	
	# coupled first-order IC and RHS for hyperbolic problems
	# initial condition (IC)
	i((u, v),) = (
		(-u, bc(u; static = true)),
		 -v
	)
	# ODE RHS
	f((u, v),) = (
		(divergence(w(v)) / p.h, bc(u)),
			    _grad(u)  / p.h
	)
	
	# second-order IC and RHS for parabolic problems
	i(u::AbstractObject) = i((u,       0*v),      )[1]
	f(u::AbstractObject) = f((u, f((u, 0*v),)[2]),)[1]

	if     imex == :explicit
		f_ex = f
		f_im = w -> 0w
	elseif imex == :implicit
		f_ex = w -> 0w
		f_im = f
	else
		error("argument 'imex' is neither of (:implicit, :explicit)")
	end

	if     bolicity == :parabolic
		intg =       tr_bdf2(f_ex, f_im,  u    )
	elseif bolicity == :hyperbolic
		intg = tr_bdf2_schur(f_ex, f_im, (u, v))
	else
		error("argument 'bolicity' is neither of (:parabolic, :hyperbolic).")
	end
	
	init!(intg, i; newton_maxit = 30, newton_rtol = 1e-6)
	
	display_2d(ax, intg)

	# Spectral radii of the composite (2nd-order) linearized rhs operator
	(ρ_ex, ρ_im) = map((f_ex, f_im)) do f
		abs <| lanczos!(linearize(f, u); h = intg.h[1:4], maxit = 100, λtol = 1e-3)[1]
	end
	
	# By the Schur determinant theorem (ignoring BC)
	if bolicity == :hyperbolic
		(ρ_ex, ρ_im) = sqrt.((ρ_ex, ρ_im))
	end
	Meta.@show (ρ_ex, ρ_im)
	
	# time step
	ν = 0.95 # safety factor
	(dt_ex, dt_im) = ν .* sqrt(3) ./ (ρ_ex, ρ_im)
	intg.dt[] = min(dt_ex, dt_im)
	Meta.@show (dt_ex, dt_im)
	
	while intg.t[] < duration
		println("STEP $i, implicit error: $(
			err = step!(intg; newton_maxit = 3, newton_rtol = 1e-6)
		)")
		
		display_2d(ax, intg)

		# dimensionless error measure
		η = err / (rtol * l2(u) + atol)

		# amend time step
		intg.dt[] = min(dt_ex, ν*intg.dt[]/cbrt(η))
		Meta.@show intg.dt[]
	end
end

using OffsetArrays

function brusselator(;a, b, x_0, y_0, nsteps, duration)
	t = OffsetArray(zeros(nsteps+1), 0:nsteps)
	x = OffsetArray(zeros(nsteps+1), 0:nsteps)
	y = OffsetArray(zeros(nsteps+1), 0:nsteps)
	
	X = Field((1,), ((1,),))
	Y = Field((1,), ((1,),))

	f((x, y),) = (
		a + x^2*y - b*x - x,
		  - x^2*y + b*x
	)
	
	assign!(X, x_0)
	assign!(Y, y_0)

	intg = tr_bdf2_schur(f, xy -> 0*xy, (X, Y))

	k    = 0
	t[k] = 0.
	x[k] = X.data[1]
	y[k] = Y.data[1]
	println("t = $(t[k]), x = $(x[k]), y = $(y[k])")

	while k < nsteps && t[k] < duration
		k += 1

		ρ = gershgorin!(linearize(f, (X, Y)), intg.g)
		
		ν = 0.1 # safety factor
		intg.dt[] = ν * sqrt(3) / abs(ρ[1])
		print("t = $(t[k-1]), dt = $(intg.dt[]), ρ = $ρ, ")

		step!(intg; newton_maxit = 0, newton_rtol = 0, quiet = true)

		t[k] = intg.t[]
		x[k] = X.data[1]
		y[k] = Y.data[1]
		println("x = $(x[k]), y = $(y[k])")

		display <| plot(t[0:k], [x[0:k] y[0:k]])
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
	
	# lid_driven(:scalar, :parabolic, :explicit, p, ax)
	brusselator(a = 1, b = 3, x_0 = 1, y_0 = 1, nsteps = 10000, duration = 30)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end