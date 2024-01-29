#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

function display_2d(ax, intg::tr_bdf2{Intg_})
	plt1 = heatmap(ax..., intg.y,    "y",    c = :davos)
	plt2 = heatmap(ax..., intg.e_ex, "e_ex", c = :davos)
	plt3 = heatmap(ax..., intg.e_im, "e_im", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1,3))
	display(plt)
end

function display_2d(ax, intg::tr_bdf2{Intg_SCR})
	plt1 = heatmap(ax..., intg.y[1], "y", c = :davos)
	plt2 = heatmap(ax..., intg.e_ex[1], "e_ex", c = :davos)
	plt3 = heatmap(ax..., intg.e_im[1], "e_im", c = :davos)
	plt  = plot(plt1, plt2, plt3; layout = (1,3))
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

using Traceur

function lid_driven(rank, bolicity, imex, p, ax; duration = 0.005, rtol = 1e-4, atol = 1e-4)
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
	
	# expr = (divergence(w(_grad(u))) / p.h, bc(u))
	# for _ in 1:50000
	# 	assign!(u, expr)
	# 	# assign!(v, _grad(u))
	# end
	# return

	# second-order IC and RHS for parabolic problems
	i(u::AbstractObject) = i((u,       0*v),      )[1]
	f(u::AbstractObject) = f((u, f((u, 0*v),)[2]),)[1]

	if     imex == :explicit
		f_ex = f
		f_im = nullfunc
	elseif imex == :implicit
		f_ex = nullfunc
		f_im = f
	else
		error("argument 'imex' is neither of (:implicit, :explicit)")
	end

	if     bolicity == :parabolic
		intg = tr_bdf2(Intg_,     u,     f_ex = f_ex, f_im = f_im)
	elseif bolicity == :hyperbolic
		intg = tr_bdf2(Intg_SCR, (u, v), f_ex = f_ex, f_im = f_im)
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
	
	k = 1
	while intg.t[] < duration
		
		println("STEP $k, t = $(intg.t[]), dt = $(intg.dt[]), ε_im,ex = $(
			((_, ε_ex, ε_im) = step!(intg; rtol = rtol, atol = atol, newton_maxit = 3, newton_rtol = 1e-6, quiet = (imex == :explicit)))[2:3]
		)")
		
		mod(k, 50) == 0 && display_2d(ax, intg)

		# amend time step
		dt_expl = intg.dt[] / sqrt(ε_ex)
		dt_impl = intg.dt[] / cbrt(ε_im)
		dt_grow = intg.dt[] / .95
		dt_stab = 0.95 * sqrt(3) / ρ_ex
		intg.dt[] = min(dt_grow, dt_expl, dt_impl, dt_stab)
		
		k += 1
	end
end

using OffsetArrays

function brusselator(; imex, a, b, x_0, y_0, nsteps, duration, rtol = 0., atol = 0.)
	t       = OffsetArray(zeros(nsteps+1), 0:nsteps)
	dt      = OffsetArray(zeros(nsteps  ), 1:nsteps)
	dt_expl = OffsetArray(zeros(nsteps  ), 1:nsteps)
	dt_impl = OffsetArray(zeros(nsteps  ), 1:nsteps)
	dt_stab = OffsetArray(zeros(nsteps  ), 1:nsteps)
	x       = OffsetArray(zeros(nsteps+1), 0:nsteps)
	y       = OffsetArray(zeros(nsteps+1), 0:nsteps)
	ε_expl  = OffsetArray(zeros(nsteps  ), 1:nsteps)
	ε_impl  = OffsetArray(zeros(nsteps  ), 1:nsteps)
	
	X = Field((1,), ((1,),))
	Y = Field((1,), ((1,),))

	assign!(X, x_0)
	assign!(Y, y_0)

	f((x, y),) = (
		a + x^2*y - b*x - x,
		  - x^2*y + b*x
	)
	
	if     imex == :explicit
		f_ex = f
		f_im = w -> 0w
	elseif imex == :implicit
		f_ex = w -> 0w
		f_im = f
	else
		error("argument 'imex' is neither of (:implicit, :explicit)")
	end

	intg = tr_bdf2_schur(f_ex, f_im, (X, Y))

	k    = 0
	t[k] = 0.
	x[k] = X.data[1]
	y[k] = Y.data[1]
	println("t = $(t[k]), x = $(x[k]), y = $(y[k])")

	while k < nsteps && t[k] < duration
		k += 1

		(ρ_ex, ρ_im) = map((f_ex, f_im)) do f
			abs <| gershgorin!(linearize(f, (X, Y)), intg.g)[1]
		end
		
		if k == 1
			dt_impl[k] = Inf
			dt_expl[k] = Inf
			dt_grow    = Inf
			dt_stab[k] = 0.05 * sqrt(3) / max(ρ_ex, ρ_im)
		else
			dt_expl[k] = dt[k-1] / sqrt(ε_expl[k-1])
			dt_impl[k] = dt[k-1] / cbrt(ε_impl[k-1])
			dt_grow    = dt[k-1] / .95
			dt_stab[k] = 0.95 * sqrt(3) / max(ρ_ex)
		end
		intg.dt[] = dt[k] = min(dt_grow, dt_expl[k], dt_impl[k], dt_stab[k])
		print("t = $(t[k-1]), dt = $(dt[k]), ")

		# step computes dimensionless error
		(ε_expl[k], ε_impl[k]) = step!(intg; rtol = rtol, atol = atol, newton_maxit = 30, newton_rtol = 1e-6, quiet = true)
		print("ε = ($(ε_expl[k]), $(ε_impl[k])), ")
		
		t[k] = intg.t[]
		x[k] = X.data[1]
		y[k] = Y.data[1]
		# println("x = $(x[k]), y = $(y[k])")
		println("")

		display <| plot(
			plot(t[0:k  ], [x[0:k]  y[0:k]] ; label=["x (-)" "y (-)"]),
			plot(t[1:k  ], [log10.( ε_expl[1:k]) log10.( ε_impl[1:k])]; label=["log₁₀ ε_ex (-)" "log₁₀ ε_im (-)"]),
			plot(t[0:k-1], [log10.(dt_expl[1:k]) log10.(dt_impl[1:k]) log10.(dt_stab[1:k]) log10.(dt[1:k])]; label=["log₁₀ dt_expl (-)" "log₁₀ dt_impl (-)" "log₁₀ dt_stab (-)" "log₁₀ dt (-)"]),
			; layout = (3,1)
		)
	end
end

function brusselator_diffusion(p, ax; a, b, x_0, y_0, Dx, Dy, nsteps, duration, rtol = 0., atol = 0.)
	t       = OffsetArray(zeros(nsteps+1), 0:nsteps)
	dt      = OffsetArray(zeros(nsteps  ), 1:nsteps)
	dt_expl = OffsetArray(zeros(nsteps  ), 1:nsteps)
	dt_impl = OffsetArray(zeros(nsteps  ), 1:nsteps)
	dt_stab = OffsetArray(zeros(nsteps  ), 1:nsteps)
	ε_expl  = OffsetArray(zeros(nsteps  ), 1:nsteps)
	ε_impl  = OffsetArray(zeros(nsteps  ), 1:nsteps)
	
	x = Field(p.n, div_stags)
	y = Field(p.n, div_stags)

	σ = 2

	assign!(x, fieldgen((i, j) -> max(0, min(4.5, x_0 + (rand()-0.5)*2*σ))))
	assign!(y, fieldgen((i, j) -> max(0, min(4.5, y_0 + (rand()-0.5)*2*σ))))

	f_ex((x, y),) = (
		(divergence(Dx*grad(x)) + a + x^2*y - b*x - x, zero_bc(0)(x)),
		(divergence(Dy*grad(y))     - x^2*y + b*x,     zero_bc(0)(y))
	)

	f_im = xy -> 0*xy
	# f_im((x, y),) = (
	# 	(divergence(Dx*grad(x)), zero_bc(0)(x)),
	# 	(divergence(Dy*grad(y)), zero_bc(0)(y))
	# )
	
	intg = tr_bdf2(f_ex, f_im, (x, y))

	k    = 0
	t[k] = 0.

	dt_im = 0

	while k < nsteps && t[k] < duration
		k += 1

		ρ_ex = abs <| gershgorin!(linearize(f_ex, (x, y)), intg.h[1:2])[1]
		
		if k == 1
			dt_impl[k] = Inf
			dt_expl[k] = Inf
			dt_grow    = Inf
			dt_stab[k] = 0.05 * sqrt(3) / ρ_ex
		else
			dt_expl[k] = dt[k-1] / sqrt(ε_expl[k-1])
			dt_impl[k] = dt[k-1] / cbrt(ε_impl[k-1])
			dt_grow    = dt[k-1] / .95
			dt_stab[k] = 0.95 * sqrt(3) / ρ_ex
		end
		intg.dt[] = dt[k] = min(dt_grow, dt_expl[k], dt_impl[k], dt_stab[k])
		println("t = $(t[k-1]), dt = $(dt[k]), ")

		# step, computes dimensionless error measure
		(t[k], ε_expl[k], ε_impl[k]) = step!(intg; rtol = rtol, atol = atol, newton_maxit = 30, newton_rtol = 1e-5, quiet = true)
		
		display <| plot(
			plot(heatmap(ax..., intg.y[1],   "x", c = :davos), heatmap(ax..., intg.y[2],   "y", c = :davos); layout = (1,2)),
			# plot(heatmap(ax..., intg.e[1], "e_x", c = :davos), heatmap(ax..., intg.e[2], "e_y", c = :davos); layout = (1,2)),
			plot(t[1:k  ], [log10.( ε_expl[1:k]) log10.( ε_impl[1:k])]; label=["log₁₀ ε_ex (-)" "log₁₀ ε_im (-)"]),
			plot(t[0:k-1], [log10.(dt_expl[1:k]) log10.(dt_impl[1:k]) log10.(dt_stab[1:k]) log10.(dt[1:k])]; label=["log₁₀ dt_expl (-)" "log₁₀ dt_impl (-)" "log₁₀ dt_stab (-)" "log₁₀ dt (-)"]),
			; layout = (3,1)
		)
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
	
	lid_driven(:vector, :parabolic, :explicit, p, ax)
	# brusselator_diffusion(p, ax; a = 1, b = 3, x_0 = 1, y_0 = 1, Dx = 0.2, Dy = 0.02, nsteps = 30000, duration = Inf, rtol = 1e-3)
	# brusselator(; imex = :explicit, a = 1, b = 3, x_0 = 1.5, y_0 = 1.5, nsteps = 30000, duration = 30, rtol = 1e-4)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end