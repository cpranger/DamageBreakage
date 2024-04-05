#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

function display_2d(k, ax, intg, (v, e, α), (j1, j2, j3), t, ε_ex, ε_im)
    plt1 = plot(
        heatmap(ax..., v.x, "v_x", c = :broc,  clims=((-1,1) .* 1.1 .* absmax(v.x))),
        heatmap(ax..., v.y, "v_y", c = :broc,  clims=((-1,1) .* 1.1 .* absmax(v.y))),
        heatmap(ax..., α,   "α",   c = :davos, clims=(( 0,1) .* 1.1 .*    max(α  ))),
        ; layout = (1,3)
    )
    # plt2 = plot(
    #     heatmap(ax..., e.xx, "e_xx", c = :broc, clims=((-1,1) .* 1.1 .* absmax(e.xx))),
    #     heatmap(ax..., e.xy, "e_yx", c = :broc, clims=((-1,1) .* 1.1 .* absmax(e.xy))),
    #     heatmap(ax..., e.yy, "e_yy", c = :broc, clims=((-1,1) .* 1.1 .* absmax(e.yy))),
    #     ; layout = (1,3)
    # )
	plt2 = plot(
        heatmap(ax..., j1, "J_1", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j1))=#),
        heatmap(ax..., j2, "J_2", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        # heatmap(ax..., j3, "J_3", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j3))=#),
        ; layout = (1,2)
    )
    # plt2 = plot(
    #     # heatmap(ax..., intg.w_1.v.x, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j1))=#),
    #     # heatmap(ax..., intg.w_1.v.y, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
    #     heatmap(ax..., intg.w_1.e.xx, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
    #     heatmap(ax..., intg.w_1.e.xy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
    #     heatmap(ax..., intg.w_1.e.yy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
    #     ; layout = (1,3)
    # )
    plt3 = plot(
        # heatmap(ax..., intg.w_2.v.x, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j1))=#),
        # heatmap(ax..., intg.w_2.v.y, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.w_2.e.xx, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.w_2.e.xy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.w_2.e.yy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        ; layout = (1,3)
    )
    plt4 = plot(
        # heatmap(ax..., intg.w_3.v.x, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j1))=#),
        # heatmap(ax..., intg.w_3.v.y, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.w_3.e.xx, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.w_3.e.xy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.w_3.e.yy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        ; layout = (1,3)
    )
    plt5 = plot(
        # heatmap(ax..., intg.y.v.x, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j1))=#),
        # heatmap(ax..., intg.y.v.y, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.y.e.xx, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.y.e.xy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.y.e.yy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        ; layout = (1,3)
    )
    plt6 = plot(
        # heatmap(ax..., intg.y.v.x, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j1))=#),
        # heatmap(ax..., intg.y.v.y, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.e_ex.e.xx, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.e_ex.e.xy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.e_ex.e.yy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        ; layout = (1,3)
    )
    plt7 = plot(
        # heatmap(ax..., intg.y.v.x, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j1))=#),
        # heatmap(ax..., intg.y.v.y, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.e_im.e.xx, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.e_im.e.xy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        heatmap(ax..., intg.e_im.e.yy, "", c = :davos#=, clims=((-1,1) .* 1.1 .* absmax(j2))=#),
        ; layout = (1,3)
    )
	plt3 = plot(
        t[1:k],
        reshape([
            map(e -> log10(e.v), ε_ex[1:k])
            map(e -> log10(e.v), ε_im[1:k])
            map(e -> log10(e.e), ε_ex[1:k])
            map(e -> log10(e.e), ε_im[1:k])
            map(e -> log10(e.α), ε_ex[1:k])
            map(e -> log10(e.α), ε_im[1:k])
            log10.(diff(t[0:k]))
        ], k, :);
        label = reshape([
            "log₁₀ ε_ex.v (-)",
            "log₁₀ ε_im.v (-)",
            "log₁₀ ε_ex.e (-)",
            "log₁₀ ε_im.e (-)",
            "log₁₀ ε_ex.α (-)",
            "log₁₀ ε_im.α (-)",
            "log₁₀ dt (-)"
        ], 1, :)
    )
	plt  = plot(plt1, plt2, plt3; layout = (3,1))
    # plt  = plot(plt2, plt3, plt4, plt5; layout = (4,1))
    # plt  = plot(plt6, plt7; layout = (2,1))
	display(plt)
end

function damage_isotropic(p, ax; nsteps, duration, rtol, atol = 0)
	t    = OffsetArray(zeros(nsteps+1), 0:nsteps)
	ε_ex = OffsetArray(fill((v = 0., e = 0., α = 0.), nsteps), 1:nsteps)
    ε_im = OffsetArray(fill((v = 0., e = 0., α = 0.), nsteps), 1:nsteps)
	
	v0 = Tensor(p.n, motion_stags)
	e0 = Tensor(p.n, Symmetric, strain_stags)
	α0 =  Field(p.n, state_stags)

    j1 = Field(p.n, state_stags)
    j2 = Field(p.n, state_stags)
    j3 = Field(p.n, state_stags)
    
    x  = (i, j) -> i * p.h - 0.5
    y  = (i, j) -> j * p.h - 0.5

    bump(R) = (i, j) -> let r = sqrt(x(i, j)^2 + y(i, j)^2)
        r < R ? cos(.5π*r/R)^2 : 0
    end

	δ = Tensor(Identity)
    I = Tensor(Ones)
	
    assign!(v0, (; x = fieldgen(y), y = fieldgen(x)))
    assign!(e0,   I - 9*δ)
	assign!(α0,  .1 * fieldgen(bump(.25)))
    
    # assign!(j1, J1(e0))
    # assign!(j2, J2(e0))
    # # assign!(j3, J3(e0))
    # display_2d(1, ax, (v0, e0, α0), (j1, j2, j3), t, ε_ex, ε_im)

    # Meta.@show(j2)
    # return
    
    bc = v -> (;
		x = (
            BC(-2, 0),
            BC(+2, 0),
            BC(-1, 0),
            BC(+1, 0)
        ),
		y = (
            BC(-2, 0),
            BC(+2, 0),
            BC(-1, 0),
            BC(+1, 0)
        )
	)
    
    diag_grad2(α) = Tensor((
		xx = A(D(α, :x), :y)*A(D(α, :x), :y),
		xy = A(D(α, :y), :x)*A(D(α, :x), :y),
		yy = A(D(α, :y), :x)*A(D(α, :y), :x)
	), Symmetric)

    # factor 2 instead of 4) takes into account a diagonal grid factor sqrt(1/2).
    diag_lap(α) = 2*A(A(α, :x), :y) - 2*α

    λ(α) = p.λ_0
    μ(α) = p.μ_0 - α * p.μ_r
    γ(α) = α * p.γ_r

	s_elast(  e, α) = (λ(α)*J1(e) - γ(α)*sqrt(J2(e))) * δ + (2*μ(α) - γ(α)*J1(e)/sqrt(J2(e))) * e
	s_struc(     α) = #=.....=# diag_grad2(α) / p.h^2
	
	f_v_ex(v, e, α) = (divergence(s_elast(e, α ) - s_elast(e, α0) + s_struc(α)) / p.h, bc(v))
	f_v_im(v, e   ) = (divergence(s_elast(e, α0)) / p.h,                               bc(v))
	
	f_e_ex(       ) = 0*δ
	f_e_im(v,     ) = symgrad(v) / p.h
	
	f_α_ex(   e   ) = p.C * p.μ_r * J2(e) + p.C * p.γ_r * J1(e) * sqrt(J2(e))
	f_α_im(      α) = p.C * p.κ * diag_lap(α) / p.h^2
	
	f_ex  =  ((v, e, α),) -> (
        f_v_ex(v, e, α),
        f_e_ex(       ),
        f_α_ex(   e   )
    )

    f_im  =  ((v, e, α),) -> (
        f_v_im(v, e   ),
        f_e_im(v      ),
        f_α_im(      α)
    )

    evo = tr_bdf2(Intg_DR, (v = v0, e = e0, α = α0), f_ex = f_ex, f_im = f_im)
    
    for k in 1:nsteps
		(ε_ex[k], ε_im[k]) = step!(evo; atol = atol, rtol = rtol, newton_maxit = 30, newton_rtol = 1e-5, quiet = false)
        t[k] = evo.t[]
        assign!(j1, J1(e0))
        assign!(j2, J2(e0))
        # assign!(j3, J3(e0))
        display_2d(k, ax, evo, (v0, e0, α0), (j1, j2, j3), t, ε_ex, ε_im)
        t[k] > duration && break 
	end
end

function parameters(; nb)
	n    =  (100, 100) .* BLOCK_SIZE   # mesh resolution
	h    =  1 / (n[1] - 2)
	l    =  (n .- 2) .* h            # physical domain size
	
	r_0  =  1
	μ_0  =  1
	λ_0  =  μ_0                      # assuming Poisson's ratio ν = 1/4
	k_0  =  λ_0 + (2/3)*μ_0          # Bulk modulus

    μ_r = .9μ_0 # damage modulus (?)
    γ_r = .1    # damage modulus (?)

	C   = 1  # damage rate constant
	κ   = 1  # diffusivity

	
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
	damage_isotropic(p, ax; nsteps = 2, duration = Inf, rtol = 1e-3)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end