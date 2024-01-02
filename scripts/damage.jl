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

function damage_isotropic(p, ax; nsteps, duration, rtol)
	t  = OffsetArray(zeros(nsteps+1), 0:nsteps)
	η  = OffsetArray(zeros(nsteps  ), 1:nsteps)
	
	v0 = Vector(p.n, motion_stags)
	e0 = Tensor(p.n, strain_stags)
	α0 =  Field(p.n, state_stags)
    
	assign!(v0, _____)
	assign!(e0, _____)
	assign!(α0, _____)

    δ = Tensor(Identity)
    Lap(x) = divergence(grad(x))
    
    a(e, α) =   λ(α) - γ(α)/J1(e)*sqrt(J2(e))
    b(e, α) = 2*μ(α) - γ(α)*J1(e)/sqrt(J2(e))
    c(   α) = 1
    
    g() = p.C * p.μ_r
    h() = p.C * p.γ_r
    i() = p.C * p.κ

    s_elast(  e, α) = a(e, α) * J1(e) * δ + b(e, α) * e
    s_struc(     α) = c(   α) * grad2(α) / p.h^2
    
    f_v_ex(   e, α) = divergence(s_elast(e, α ) - s_elast(e, α0) + s_struc(α)) / p.h
    f_v_im(   e   ) = divergence(s_elast(e, α0)) / p.h
    
	f_e_ex(       ) = 0*δ
    f_e_im(v,     ) = symgrad(v) / p.h
    
    f_α_ex(   e   ) = g() * J2(e) + h() * J1(e) * sqrt(J2(e))
    f_α_im(      α) = i() * Lap(α) / p.h^2

    f_ex   =  (_, e, α) -> (
        f_v_ex(   e, α),
        f_e_ex(       ),
        f_α_ex(   e   )
    )

    f_im   =  (v, e, α) -> (
        f_v_im(   e   ),
        f_e_im(v      ),
        f_α_im(      α)
    )

    evo = tr_bdf2_dr(f_ex, f_im, (v0, e0, α0), Ref(1.))
    
    k    = 0
	t[k] = 0.
    
	while k < nsteps && t[k] < duration
		k += 1

		# ρ_ex = abs <| gershgorin!(
        #     linearize(
        #         ((e , α ),) -> (f_e_ex(       ), f_α_ex(   e   )), # TODO: check eigenvalues analytically
        #          (e0, α0)
        #     ),
        #     ___ # aux vecs # TODO
        # )[1]
        
		# step computes dimensionless error measure
		(t[k], η_ex[k], η_im[k]) = step!(evo; atol = atol, rtol = rtol, newton_maxit = 30, newton_rtol = 1e-5, quiet = false)
		
		display <| plot(
			plot(heatmap(ax..., evo.y.v,  "v", c = :davos), heatmap(ax..., evo.y.α, "α", c = :davos); layout = (1,2)),
			plot(heatmap(ax..., evo.α[1], "e_x", c = :davos), heatmap(ax..., evo.e[2], "e_y", c = :davos); layout = (1,2)),
			plot(t[1:k], [log10.(η[1:k]) log10.(diff(t))]; label = ["log₁₀ ε (-)" "log₁₀ dt (-)"]),
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
	
	# lid_driven(:scalar, :parabolic, :explicit, p, ax)
	damage_isotropic(p, ax; nsteps = 30000, duration = Inf, rtol = 1e-3)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end