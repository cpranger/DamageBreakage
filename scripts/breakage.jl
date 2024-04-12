#!/usr/bin/env -S julia --project

include("./header.jl")

using SpecialFunctions
using LambertW

using StaggeredKernels.Antiline

function display_1d(k, p, x, (v, e, β), t, ε)
    plts = (
        plot(x.data[2, 1:end],  v.y.data[1, 1:end]; label="v (m/s)"),
        plot(x.data[1, 1:end-1], e.xy.data[1, 1:end-1]; label="e (-)"),
        plot(x.data[1, 1:end-1],    β.data[1, 1:end-1]; label="γ (-)"),
        plot(t[1:k], [log10.(diff(t[0:k])) log10.(ε[1:k])]; label=["log₁₀ dt (-)" "log₁₀ ε (-)"])
    )
    display <| plot(plts...; layout=(4,1))
end

# Based on Pranger et al., 2022, cycle model
function breakage_1d(p; nsteps, duration, rtol, atol = 0)
	t = OffsetArray(zeros(nsteps+1), 0:nsteps)
    ε = OffsetArray(zeros(nsteps), 1:nsteps)
	# ε = OffsetArray(fill((v = 0., e = 0., β = 0.), nsteps), 1:nsteps)
	
	x  =  Field(p.n, ((0,), (1,)))
    v0 = Tensor(p.n, motion_stags)
	e0 = Tensor(p.n, Symmetric, strain_stags)
	β0 =  Field(p.n, state_stags)
    γ  =  Field(p.n, state_stags)

    assign!(x,  fieldgen(i -> i*p.h)) # x-axis
	assign!(v0, p.v_init(x))
    assign!(e0, p.e_init(x))
	assign!(β0, p.β_init(x))
    assign!(γ,  p.γ(β0))

    display_1d(0, p, x, (v0, e0, γ), t, ε)
    
    f_ex  =  ((v, e, β),) -> (
        p.f_v_ex(v      ),
        p.f_e_ex(      β),
        p.f_β_ex(v, e, β)
    )

    f_im  =  ((v, e, β),) -> (
        p.f_v_im(v, e   ),
        p.f_e_im(v      ),
        p.f_β_im(      β)
    )

    evo = tr_bdf2(Intg_DR, (v = v0, e = e0, α = β0), f_ex = f_ex, f_im = f_im)
    
    for k in 1:nsteps
		ε[k] = step!(evo; atol = atol, rtol = rtol, newton_maxit = 30, newton_rtol = 1e-5, quiet = true)
        t[k] = evo.t[]
        assign!(γ,  p.γ(β0))
        display_1d(k, p, x, (v0, e0, γ), t, ε)
        t[k] > duration && break 
	end
end

function parameters(#=; nb=#)
	n    =  (100,)                  # mesh resolution
	
    a    =  0.02                    # rate-and-state dimensionless parameter
	b    =  0.03                    # rate-and-state dimensionless parameter

	σ_n  =  1e6                     # normal stress, Pa
	V_0  =  1e-6                    # reference velocity, m/s
	V_p  =  1e-9                    # driving velocity, m/s
	r_0  =  1e3                     # density, kg/m^3
	μ_0  =  1e10                    # shear modulus, Pa
	λ_0  =  μ_0                     # Lamé parameter, assuming Poisson's ratio ν = 1/4
	k_0  =  λ_0 + (2/3)*μ_0         # bulk modulus, Pa
    c_s  =  sqrt(μ_0 / r_0)         # solenoidal  wave speed, m/s
	c_p  =  sqrt((λ_0+2*μ_0)/r_0)   # compressive wave speed, m/s
	
    d_c  =  0.01                    # slip weakening distance, m
	H    =  50000                   # boundary distance, m
	w    =  0.1                     # model domain size, m
    h    =  w / (n[1] - 2)          # cell size, m
	
	λ    =  18 * h                  # diffusion length scale, m
	
	η    =  1/2 * μ_0/σ_n * V_0/c_s # radiation damping viscosity, -
    
	c_2  =  1 + (sqrt(2) - 1) * a/b # dimensionless coefficient 1
	c_1  =  2/(1 + a/b) * c_2       # dimensionless coefficient 2
	c_3  =  c_1                     # dimensionless coefficient 3
    c_4  =  1/(π*(1-a/b))           # dimensionless coefficient 4
    
	γ_0     = (V_0/(sqrt(π)*2*λ)) * gamma(1 + 1/(1 - a/b)) / gamma(1/2 + 1/(1 - a/b)) # reference strain rate, 1/s
	ψ_0     =  V_0/d_c              # reference state
	
    γ(β)    =  a/(c_3*η) * lambertw((c_3*η)/a * exp(β/a))
	ψ(β, f) = γ(β)^(a/b) * exp(-(1/b)*(f - c_3*η*γ(β)))

	v_init(x) = V_p / H * x
    e_init(x) = 0.16 * σ_n/μ_0
    β_init_(x) = 1e-4 + If(abs(c_4*x/λ) < 1, 2*c_4*cos(1/2*π*c_4*x/λ)^2, 0)
    β_init(x)  = c_3*η*β_init_(x) + a*log(β_init_(x))
    
	bc_v = v -> (;
        y = (
            BC(-1, 0),
            BC(+1, 0)
        )
	)
    
    bc_β = β -> (
		BC(-1, FD(β, :x)),
        BC(+1, BD(β, :x))
	)
    
    f_v_ex(v      ) = 0*v
	f_v_im(v, e   ) = c_s^2 * (divergence(e), bc_v(v)) / h
	
	f_e_ex(      β) = - γ_0 * γ(β) * (; xy = 1)
	f_e_im(v,     ) = symgrad(v) / h
	
    f_β_ex(v, e, β) = - b * c_2 * ψ_0 * ψ(β, μ_0*J2(e)/σ_n) +
                      + b * c_1 * ψ_0 * γ(β) +
                      + μ_0/σ_n * J2(symgrad(v))/h
	f_β_im(      β) =   b * c_1 * ψ_0 * (1-a/b)^2 * λ^2 * (divergence(grad(β)), bc_β(β)) / h^2

    # collect all variables local to this function:
	vars = Base.@locals
	
	# and return as a named tuple
	return NamedTuple{Tuple(keys(vars))}(values(vars))
end

function main()
	s = ArgParseSettings()
	
	@add_arg_table s begin
    	# "--nb"
		# 	help = "comma-separated list of the number of blocks in each dimension"
		# 	required = true
		# 	arg_type = (NTuple{N, Int} where N)
    end
	
	p = parameters(; parse_args(s)...);# pprintln(p)
	
	breakage_1d(p; nsteps = 100000, duration = Inf, rtol = 1e-5)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end