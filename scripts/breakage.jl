#!/usr/bin/env -S julia --project

include("./header.jl")

using SpecialFunctions
using LambertW

using StaggeredKernels.Antiline

function display_1d(k, p, x, (u, v, e, β), t, dt_, ε)
    plts = (
        plot(x.data[2, 1:end],  u.y.data[1, 1:end]; label="u (m)"),
        plot(x.data[2, 1:end],  v.y.data[1, 1:end]; label="v (m/s)"),
        plot(x.data[1, 1:end-1], e.xy.data[1, 1:end-1]; label="e (-)"),
        plot(x.data[1, 1:end-1],    β.data[1, 1:end-1]; label="γ (-)"),
        # plot(t[1:k] .- t[1], [log10.(diff(t[0:k])) log10.(dt_[1:k])]; label=["log₁₀ dt (-)" "log₁₀ dt* (-)"]),
        plot(t[max(k-50, 1):k] .- t[max(k-50, 1)], [log10.(diff(t[max(k-51, 0):k])) log10.(dt_[max(k-50, 1):k])]; label=["log₁₀ dt (-)" "log₁₀ dt* (-)"]),
        # plot(t[1:k] .- t[1], [log10.(ε.e for ε in ε[1:k]) log10.(ε.α for ε in ε[1:k])]; label=["log₁₀ ε.e (-)" "log₁₀ ε.β (-)"]),
        plot(t[max(k-50, 1):k] .- t[max(k-50, 1)], [log10.(ε.u for ε in ε[max(k-50, 1):k]) log10.(ε.v for ε in ε[max(k-50, 1):k]) log10.(ε.e for ε in ε[max(k-50, 1):k]) log10.(ε.α for ε in ε[max(k-50, 1):k])]; label=["log₁₀ ε.u (-)" "log₁₀ ε.v (-)" "log₁₀ ε.e (-)" "log₁₀ ε.β (-)"])
    )
    display <| plot(plts...; layout=(length(plts),1))
end

# Based on Pranger et al., 2022, cycle model
function breakage_1d(p; nsteps, duration, rtol, atol = 0)
	t   = OffsetArray(zeros(nsteps+1), 0:nsteps)
    dt_ = OffsetArray(zeros(nsteps  ), 1:nsteps)
    ε   = OffsetArray(fill((u = 0., v = 0., e = 0., α = 0.), nsteps), 1:nsteps)
    
	x  =  Field(p.n, ((0,), (1,)))
    u0 = Tensor(p.n, motion_stags)
	v0 = Tensor(p.n, motion_stags)
	e0 = Tensor(p.n, Symmetric, strain_stags)
	β0 =  Field(p.n, state_stags)
    γ  =  Field(p.n, state_stags)

    assign!(x,  p.x_init())
    # assign!(u0, 0) # p.u_init(x)
    # assign!(v0, p.v_init(x))
    assign!(β0, p.β_init(x))
    assign!(γ,  p.γ_0 * p.γ(β0))

    assign!(v0.y, (B(B(v0.y, :x), :x) + B(γ, :x) * p.dx, (BC(-1, 0),)))
    assign!(v0.y, (1-x/p.d_d) * v0.y + (x/p.d_d) * p.v_d)

    assign!(e0, p.e_0)

    assign!(u0.y, p.e_0 * x)

    p.t[] = p.d_d * p.e_0 / p.v_d

    Meta.@show p.t[] / pi / 1e7

    display_1d(0, p, x, (u0, v0, e0, γ), t, dt_, ε)
    
    f_ex  =  (t, (u, v, e, β),) -> (
        p.f_u_ex(t,    v      ),
        p.f_v_ex(t,    v      ),
        p.f_e_ex(t,          β),
        p.f_β_ex(t,    v, e, β)
    )

    f_im  =  (t, (u, v, e, β),) -> (
        p.f_u_im(t, u         ),
        p.f_v_im(t, u,    e   ),# - 2*p.c_σ/p.dx * v,
        p.f_e_im(t,    v      ),
        p.f_β_im(t,          β)
    )

    evo = tr_bdf2(Intg_DR, (u = u0, v = v0, e = e0, α = β0), f_ex = f_ex, f_im = f_im, t = p.t, dt = p.dt, dt_ = p.dt_)
    
    v_atol = 0.
    v_rtol = rtol

    for k in 1:nsteps
		ε[k] = step!(
            evo;
            atol = (; u = atol, v = v_atol, e = atol, α = atol),
            rtol = (; u = rtol, v = v_rtol, e = rtol, α = rtol),
            newton_maxit = 30,
            newton_rtol = 1e-4,
            cg_maxit = Int <| round <| 10*sqrt(sum(p.n.^2)),
            cg_rtol = 1e-5,
            growth = 1.1,
            safety = 1
        )
        t[k]   = evo.t[]
        dt_[k] = evo.dt_[]
        assign!(γ,  p.γ(β0))
        display_1d(k, p, x, (u0, v0, e0, γ), t, dt_, ε)
        v_rtol = min(1.1*v_rtol, 1e6)
        t[k] > duration && break
	end
end

function parameters()
	n    =  (200,)                  # mesh resolution
	
    a    =  0.02                    # rate-and-state dimensionless parameter
	b    =  0.03                    # rate-and-state dimensionless parameter
    
	σ_n  =  1e6                     # normal stress, Pa
	v_0  =  1e-6                    # reference velocity, m/s
    v_d  =  1e-9                    # driving   velocity, m/s
    r_0  =  1e3                     # density, kg/m^3
	μ_0  =  1e10                    # shear modulus, Pa
	λ_0  =  μ_0                     # Lamé parameter, assuming Poisson's ratio ν = 1/4
	k_0  =  λ_0 + (2/3)*μ_0         # bulk modulus, Pa
    c_s  =  sqrt(μ_0 / r_0)         # solenoidal  wave speed, m/s
	c_p  =  sqrt((λ_0+2*μ_0)/r_0)   # compressive wave speed, m/s
	
    d_d  =  5e4                     # driving distance, m
	d_c  =  1e-2                    # slip weakening distance, m
    dx   =  1e-4                    # cell size
    w    =  (n[1]-2) * dx           # domain size
    λ    =  18 * dx                 # diffusion length scale, m
	
	η    =  1/2 * μ_0/σ_n * v_0/c_s # radiation damping viscosity, -
    
	c_2  =  1 + (sqrt(2) - 1) * a/b # dimensionless coefficient 2
	c_1  =  2/(1 + a/b) * c_2       # dimensionless coefficient 1
	c_3  =  c_1                     # dimensionless coefficient 3
    c_4  =  1/(π*(1-a/b))           # dimensionless coefficient 4
    
	γ_0  = (v_0/(sqrt(π)*2*λ)) * gamma(1 + 1/(1 - a/b)) / gamma(1/2 + 1/(1 - a/b)) # reference strain rate, 1/s
	ψ_0  =  v_0/d_c              # reference state
	
    e_0  =  0.16 * σ_n/μ_0
    # e_d  =  v_d / d_d

    t    = Ref(BigFloat(0.))
    dt   = Ref(1e-6)
    dt_  = Ref(0.)

    γ(β)    =  a/(c_3*η) * lambertw((c_3*η)/a * exp(β/a))
	ψ(β, f) = γ(β)^(a/b) * exp(-(1/b)*(f - c_3*η*γ(β)))

	x_init()  = fieldgen(i -> dx*i)
    u_init(x) = -x * e_0
    v_init(x) = v_d * x/d_d
    e_init(x) = e_0
    β_ini_(x) = 1e-4 + If(abs(c_4*x/λ) < 1, 2*c_4*cos(1/2*π*c_4*x/λ)^2, 0)
    β_init(x) = c_3*η*β_ini_(x) + a*log(β_ini_(x))
    
	c_σ = 1e3 * c_s
    
    lap(β) = (
        divergence(grad(β)),
        (
            BC(-1, +2*FD(β, :x)),
            BC(+1, -2*BD(β, :x))
        )
    )

    div(e, u, t) = (
        divergence(e),
        (;
            y = (
                BC(-1, 0),
                BC(+1, (v_d*t - u.y)/(d_d - w) - B(e.xy, :x))
            )
        )
    )

    f_u_ex(t,    v      ) =   v
    f_u_im(t, u         ) = 0*u
    
    f_v_ex(t,    v      ) = 0*v
    f_v_im(t, u,    e   ) =   c_σ^2 * div(e, u, t) / dx
	
	f_e_ex(t,          β) = - γ_0 * γ(β) * (; xy = 1)
	f_e_im(t,    v,     ) =   symgrad(v) / dx
	
    f_β_ex(t,    v, e, β) = - b * c_2 * ψ_0 * ψ(β, μ_0*J2(e)/σ_n) +
                            + μ_0/σ_n * J2(symgrad(v) / dx)
	f_β_im(t,          β) = + b * c_1 * ψ_0 * γ(β) +
                            + b * c_1 * ψ_0 * (1-a/b)^2 * λ^2 * lap(γ(β)) / dx^2

    # collect all variables local to this function:
	vars = Base.@locals
	
	# and return as a named tuple
	return NamedTuple{Tuple(keys(vars))}(values(vars))
end

function main()
	p = parameters()#; pprintln(p)
	
    # global out = open("out.txt", "w")
    
    global algo_depth = 0 # reset just in case
    global  verbosity = 0

	breakage_1d(p; nsteps = 100000, duration = Inf, rtol = 1e-3)

    # close(out)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end