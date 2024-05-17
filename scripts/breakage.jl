#!/usr/bin/env -S julia --project

include("./header.jl")

using SpecialFunctions
using LambertW

using StaggeredKernels.Antiline

function display_1d(i, p)
    plts = (
        plot(p.x.data[2, 1:end],    p.u.y.data[1, 1:end]  ; label="u (m)"  , legend=false),
        plot(p.x.data[2, 1:end],    p.v.y.data[1, 1:end]  ; label="v (m/s)", legend=false),
        plot(p.x.data[1, 1:end-1], p.e.xy.data[1, 1:end-1]; label="e (-)"  , legend=false),
        plot(p.x.data[1, 1:end-1],    p.γ.data[1, 1:end-1]; label="γ (-)"  , legend=false),
        # plot(i.tt[0:i.k[]] .- i.tt[0], [log10.(diff(i.tt[0:i.k[]])) log10.(diff(i.ss[0:i.k[]]))]; label=["log₁₀ dt (-)" "log₁₀ dt* (-)"]),
        plot(i.tt[max(i.k[]-51, 0):i.k[]-1] .- i.tt[max(i.k[]-51, 0)], [log10.(diff(i.tt[max(i.k[]-51, 0):i.k[]])) #=log10.(diff(i.ss[max(i.k[]-51, 0):i.k[]]))=#]; label=["log₁₀ dt (-)" #="log₁₀ dt* (-)"=#], legend=false),
        # plot(i.tt[0:i.k[]] .- i.tt[0], [log10.(ε.e for ε in i.ee[1:k]) log10.(ε.α for ε in i.ee[1:i.k[]])]; label=["log₁₀ ε.e (-)" "log₁₀ ε.β (-)"]),
        plot(i.tt[max(i.k[]-51, 0):i.k[]-1] .- i.tt[max(i.k[]-51, 0)], [log10.(ε.u for ε in i.ee[max(i.k[]-50, 1):i.k[]]) log10.(ε.v for ε in i.ee[max(i.k[]-50, 1):i.k[]]) log10.(ε.e for ε in i.ee[max(i.k[]-50, 1):i.k[]]) log10.(ε.α for ε in i.ee[max(i.k[]-50, 1):i.k[]])]; label=["log₁₀ ε.u (-)" "log₁₀ ε.v (-)" "log₁₀ ε.e (-)" "log₁₀ ε.β (-)"], legend=false)
    )
    display <| plot(plts...; layout=(length(plts),1))
end

# Based on Pranger et al., 2022, cycle model
function breakage_1d(p)
	f_ex  =  (t, (u, v, e, β),) -> (
        p.f_u_ex(t,    v      ),
        p.f_v_ex(t,    v      ),
        p.f_e_ex(t,          β),
        p.f_β_ex(t,    v, e, β)
    )

    f_im  =  (t, (u, v, e, β),) -> (
        p.f_u_im(t, u         ),
        p.f_v_im(t,    v, e   ),
        p.f_e_im(t,    v      ),
        p.f_β_im(t,          β)
    )

    i = tr_bdf2(Intg_DR, (u = p.u, v = p.v, e = p.e, α = p.β), f_ex = f_ex, f_im = f_im, dt = Ref(1e-6))
    
    display_1d(i, p)
    
    v_atol = 0.
    v_rtol = p.rtol

    for k in 1:p.nsteps
		step!(
            i;
            atol = (; u = Inf, v = v_atol, e = p.atol, α = p.atol),
            rtol = (; u = Inf, v = v_rtol, e = p.rtol, α = p.rtol),
            newton_maxit = 30,
            newton_rtol = 1e-4,
            cg_maxit = Int <| round <| 10*sqrt(sum(p.n.^2)),
            cg_rtol = 1e-5,
            growth = 1.1,
            safety = .0002,
            relax  = .9
        )
        assign!(p.γ,  p.γ_(p.β))
        display_1d(i, p)
        v_rtol = min(1.05*v_rtol, 1e6)
        i.t[] > p.duration && break
	end
end

function parameters()
	nsteps = 10000
    duration = Inf
    rtol = 1e-4
    atol = 0.
    
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
    
    dx0  =  1e-4                    # first cell size
    dx_  =  Field(n, ((1,),))       # basic grid
    dx   =  Field(n, ((0,), (1,)))  # full  grid
    x_   =  Field(n, ((0,),))       # basic grid
    x    =  Field(n, ((0,), (1,)))  # full  grid

    assign!(dx_, fieldgen(i -> dx0 * 1.00^(i+0.5)))
    assign!(dx , interpolate(dx_))

    assign!(x_, (B(B(x_, :x), :x) + B(dx_, :x), (BC(-1, 0),)))
    assign!(x,  (interpolate(x_), (BC(+1, l1(dx_)),)))
    
    w = max(x)                      # domain size
    # w    =  (n[1]-2) * dx         # domain size

    λ    =  18 * dx0                # diffusion length scale, m
	
	η    =  1/2 * μ_0/σ_n * v_0/c_s # radiation damping viscosity, -
    
	c_2  =  1 + (sqrt(2) - 1) * a/b # dimensionless coefficient 2
	c_1  =  2/(1 + a/b) * c_2       # dimensionless coefficient 1
	c_3  =  c_1                     # dimensionless coefficient 3
    c_4  =  1/(π*(1-a/b))           # dimensionless coefficient 4
    
	γ_0  = (v_0/(sqrt(π)*2*λ)) * gamma(1 + 1/(1 - a/b)) / gamma(1/2 + 1/(1 - a/b)) # reference strain rate, 1/s
	ψ_0  =  v_0/d_c                 # reference state
	
    e_0  =  0.16 * σ_n/μ_0
    
    γ_(β)    =  a/(c_3*η) * lambertw((c_3*η)/a * exp(β/a))
	ψ_(β, f) = γ_(β)^(a/b) * exp(-(1/b)*(f - c_3*η*γ_(β)))

	# x  =  Field(n, ((0,), (1,)))
    u  = Tensor(n, motion_stags)
	v  = Tensor(n, motion_stags)
	e  = Tensor(n, Symmetric, strain_stags)
	β  =  Field(n, state_stags)
    γ  =  Field(n, state_stags)

    # assign!(x,     fieldgen(i -> dx*i))
    assign!(β,     1e-4 + If(abs(c_4*x/λ) < 1, 2*c_4*cos(1/2*π*c_4*x/λ)^2, 0))
    assign!(β,     c_3*η*β + a*log(β))
    assign!(γ,     γ_0 * γ_(β))
    assign!(v.y,  (B(B(v.y, :x), :x) + 2*B(γ, :x) * dx, (BC(-1, 0),)))
    assign!(e,     e_0)
    assign!(u.y,   e_0 * x)

    v_0 = max(v.y)
    c_σ = 1e0 * c_s
    
    lap(β) = (
        divergence(grad(β)/dx)/dx,
        (
            BC(-1, +2*FD(β, :x)/dx/dx),
            BC(+1, -2*BD(β, :x)/dx/dx)
        )
    )

    div(e, v) = (
        divergence(e),
        (;
            y = (
                BC(-1, 0),
                BC(+1, -(B(e.xy, :x) - e_0) - (v.y - v_0)/c_σ)
            )
        )
    )

    f_u_ex(t,    v      ) =   v
    f_u_im(t, u         ) = 0*u
    
    f_v_ex(t,    v      ) = 0*v
    f_v_im(t,    v, e   ) =   c_σ^2 * div(e, v) / dx
	
	f_e_ex(t,          β) = - γ_0 * γ_(β) * (; xy = 1)
	f_e_im(t,    v,     ) =   symgrad(v) / dx
	
    f_β_ex(t,    v, e, β) = cos(1/2*π*x/w)^2 * (- b * c_2 * ψ_0 * ψ_(β, μ_0*J2(e)/σ_n) +
                                                + μ_0/σ_n * J2(symgrad(v) / dx))
	f_β_im(t,          β) = cos(1/2*π*x/w)^2 * (+ b * c_1 * ψ_0 * γ_(β) +
                                                + b * c_1 * ψ_0 * (1-a/b)^2 * λ^2 * lap(γ_(β)))
    
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

	breakage_1d(p)

    # close(out)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end