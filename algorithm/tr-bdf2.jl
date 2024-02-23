
const nullfunc = y -> 0*y

abstract type IntgType end

struct Intg_    <: IntgType; end  # default for parabolic problems
struct Intg_SCR <: IntgType; end  # Schur Complement Reduction for hyperbolic problems
struct Intg_DR  <: IntgType; end  # Damage Rheology
struct Intg_DBR <: IntgType; end  # Damage-Breakage Rheology

struct tr_bdf2{T <: IntgType}
	y
	w_1
	w_2
	w_3
	e_im
	e_ex
	rhs
	h   #      helper vecs
	g   # more helper vecs
	lhs
	stage_1_rhs
	stage_2_rhs
	impl_err_rhs
	expl_err
	update
	imex
	t::FieldVal
	dt::FieldVal
end

tr_bdf2_data(                  y                              ) = (; y = y, w_1 = deepcopy(y), w_2 = deepcopy(y), w_3 = deepcopy(y), e_im = deepcopy(y), e_ex = deepcopy(y), rhs = deepcopy(y))
tr_bdf2_data(::Type{Intg_},    y                              ) = (; tr_bdf2_data(y)..., h = [deepcopy(y)    for _ in 1:7],                                                                            g = [])
tr_bdf2_data(::Type{Intg_SCR}, y::Tuple                       ) = (; tr_bdf2_data(y)..., h = [deepcopy(y[1]) for _ in 1:7],                                                                            g = [deepcopy(y)           for _ in 1:2])
tr_bdf2_data(::Type{Intg_DR},  y::NamedTuple{(:v, :e, :α)}    ) = (; tr_bdf2_data(values(y))..., h = (; v = [deepcopy(y.v) for _ in 1:7], α = [deepcopy(y.α) for _ in 1:7]),                                   g = [deepcopy((; e = y.e)) for _ in 1:2])
tr_bdf2_data(::Type{Intg_DBR}, y::NamedTuple{(:v, :e, :α, :β)}) = (; tr_bdf2_data(values(y))..., h = (; v = [deepcopy(y.v) for _ in 1:7], α = [deepcopy(y.α) for _ in 1:7], β = [deepcopy(y.β) for _ in 1:7]), g = [deepcopy((; e = y.e)) for _ in 1:2])

tr_bdf2_funcs(f_im, f_ex, data, dt::FieldVal) = (
	                             tr_bdf2_lhs(               f_im,                               dt),
	prepare_assignment(data.w_1, tr_bdf2_stage_1_rhs( f_ex, f_im, data.w_1,                     dt)),
	prepare_assignment(data.w_1, tr_bdf2_stage_2_rhs( f_ex, f_im, data.w_1, data.w_2,           dt)),
	prepare_assignment(data.w_1, tr_bdf2_impl_err_rhs(      f_im, data.w_1, data.w_2, data.w_3, dt)),
	prepare_assignment(data.w_1, tr_bdf2_expl_err(    f_ex,       data.w_1, data.w_2, data.w_3, dt)),
	prepare_assignment(data.w_1, tr_bdf2_update(      f_ex, f_im, data.w_1, data.w_2, data.w_3, dt))
)

imex(f_im, f_ex) = (;
	im = typeof(f_im) != typeof(nullfunc),
	ex = typeof(f_ex) != typeof(nullfunc)
)

tr_bdf2(type, y; f_ex = nullfunc, f_im = nullfunc, t = Ref(0.), dt = Ref(0.)) = tr_bdf2_(;
	type = type,
	f_ex = f_ex,
	f_im = f_im, 
	data = tr_bdf2_data(type, y),
	dt = dt,
	t = t
)
tr_bdf2_(; type, f_ex, f_im, data, t, dt) = tr_bdf2__(; 
	type = type,
	data = data,
	funcs = tr_bdf2_funcs(f_im, f_ex, data, FieldVal(dt)),
	imex = imex(f_im, f_ex),
	dt = dt,
	t = t
)
tr_bdf2__(      ; type, data, funcs, imex, t, dt) = tr_bdf2{type}(data..., funcs..., imex, FieldVal(t), FieldVal(dt))

tr_bdf2_lhs(f_im, dt) = w -> w - dt * (1/1-sqrt(2)/2) * f_im(w)

tr_bdf2_stage_1_rhs(f_ex, f_im, w_1, dt) = w_1 + dt * (
	(2/1-sqrt(2)/1) * f_ex(w_1)
  + (1/1-sqrt(2)/2) * f_im(w_1)
)

tr_bdf2_stage_2_rhs(f_ex, f_im, w_1, w_2, dt) = w_1 + dt * (
	(1/2-sqrt(2)/3) * f_ex(w_1) + (sqrt(2)/4) * f_im(w_1)
  + (1/2+sqrt(2)/3) * f_ex(w_2) + (sqrt(2)/4) * f_im(w_2)
)

tr_bdf2_impl_err_rhs(f_im, w_1, w_2, w_3, dt) = dt * (
	( 1/3-sqrt(2)/3) * f_im(w_1)
  + ( 1/3          ) * f_im(w_2)
  + (-2/3+sqrt(2)/3) * f_im(w_3)
)

tr_bdf2_expl_err(f_ex, w_1, w_2, w_3, dt) = dt * (
	( 1/2-sqrt(2)/2) * f_ex(w_1)
  + ( 1/2          ) * f_ex(w_2)
  + (-1/1+sqrt(2)/2) * f_ex(w_3)
)

tr_bdf2_update(f_ex, f_im, w_1, w_2, w_3, dt) = w_1 + dt * (
	(    sqrt(2)/4)*(f_im(w_1) + f_ex(w_1))
  + (    sqrt(2)/4)*(f_im(w_2) + f_ex(w_2))
  + (1/1-sqrt(2)/2)*(f_im(w_3) + f_ex(w_3))
)

function init!(i::tr_bdf2{Intg_}, func; newton_maxit, newton_rtol, quiet = false)
	null = 0 * i.y
	assign!(i.y, func(null))
	newtonit!(func, i.y, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
end


function init_stepsize!(i::tr_bdf2; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, ρ_ex = 0)
	dt = min(1, sqrt(3) / ρ_ex)
	t  = i.t[]
	for _ in 1:10
		i.dt[] = dt
		(η_ex, η_im) = step!(i,
			rtol = rtol,
			atol = atol,
			newton_maxit = newton_maxit,
			newton_rtol = newton_rtol,
			quiet = quiet,
			ρ_ex = 0
		)
		i.t[]  = t
		assign!(i.y, i.w_1)
		quiet || println("dt = $dt: η = $((η_ex, η_im))")
		max(max(η_ex...)-1) < rtol && max(max(η_im...)-1) < rtol && break
		dt *= min(1/sqrt(max(η_ex...)), 1/cbrt(max(η_im...)))
	end

	i.dt[] = dt
	return dt
end


function step!(i::tr_bdf2{Intg_}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, ρ_ex = 0, growth = 1.15)
	
	i.dt[] == 0. && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)

	i.dt[] = min(i.dt[], sqrt(3) / ρ_ex)

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 stage 1:")
	if i.imex.im
		assign!(i.rhs, i.stage_1_rhs)
		assign!(i.w_2, i.rhs)
		newtonit!(w -> i.rhs - i.lhs(w), i.w_2,  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_2, i.stage_1_rhs)
	end

	quiet || println("TR-BDF2 stage 2:")
    if i.imex.im
		assign!(i.rhs, i.stage_2_rhs)
		assign!(i.w_3, i.rhs)
		newtonit!(w -> i.rhs - i.lhs(w), i.w_3,  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_3, i.stage_2_rhs)
	end

	quiet || println("TR-BDF2 explicit error:")
	assign!(i.e_ex, i.expl_err)
    
	quiet || println("TR-BDF2 implicit error:")
	if i.imex.im
		assign!(i.rhs,  i.impl_err_rhs)
		assign!(i.e_im, i.rhs)
    	newtonit!(w -> i.rhs - i.lhs(w), i.e_im, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.e_im, i.impl_err_rhs)
	end

	quiet || println("TR-BDF2 update:")
	assign!(i.y, i.update)
	
	i.t[] += i.dt[]

	η_ex = l2(i.e_ex) / (rtol * l2(i.y) + atol)
	η_im = l2(i.e_im) / (rtol * l2(i.y) + atol)

	quiet || println("t = $(i.t[]), dt = $(i.dt[]), η_ex = $η_ex, η_im = $η_im")
	
	i.t[] += i.dt[]

    i.dt[] = min(1/sqrt(max(η_ex...)), 1/cbrt(max(η_im...)), growth) * i.dt[]

	return (η_ex, η_im)
end


function init!(i::tr_bdf2{Intg_SCR}, func; newton_maxit, newton_rtol, quiet = false)
	null = 0 .* i.y
	assign!(i.y, func(null))
	sc = SchurComplement(func, i.y[2])
	newtonit!(sc, i.y[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
end

function step!(i::tr_bdf2{Intg_SCR}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, ρ_ex = 0, growth = 1.15)
	
	i.dt[] == 0. && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)

	i.dt[] = min(i.dt[], sqrt(3) / ρ_ex)

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 Schur stage 1:")
	if i.imex.im
		assign!(i.rhs, i.stage_1_rhs)
		assign!(i.w_2, i.rhs)
		sc = SchurComplement(w -> i.rhs - i.lhs(w), i.w_2[2])
		newtonit!(sc, i.w_2[1],  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_2, i.stage_1_rhs)
	end

    quiet || println("TR-BDF2 Schur stage 2:")
    if i.imex.im
		assign!(i.rhs, i.stage_2_rhs)
		assign!(i.w_3, i.rhs)
		sc = SchurComplement(w -> i.rhs - i.lhs(w), i.w_3[2])
		newtonit!(sc, i.w_3[1],  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_3, i.stage_2_rhs)
	end

	quiet || println("TR-BDF2 Schur explicit error:")
	assign!(i.e_ex, i.expl_err)
	
	quiet || println("TR-BDF2 Schur implicit error:")
	if i.imex.im
		assign!(i.rhs,  i.impl_err_rhs)
		assign!(i.e_im, i.rhs)
    	sc = SchurComplement(w -> i.rhs - i.lhs(w), i.e_im[2])
		newtonit!(sc, i.e_im[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.e_im, i.impl_err_rhs)
	end

	quiet || println("TR-BDF2 Schur update:")
	assign!(i.y, i.update)
	
	i.t[] += i.dt[]

    η_ex = l2(i.e_ex) / (rtol * l2(i.y) + atol)
	η_im = l2(i.e_im) / (rtol * l2(i.y) + atol)

	quiet || println("t = $(i.t[]), dt = $(i.dt[]), η_ex = $η_ex, η_im = $η_im")
	
	i.dt[] = min(1/sqrt(max(η_ex...)), 1/cbrt(max(η_im...)), growth) * i.dt[]

	return (η_ex, η_im)
end

function step!(i::tr_bdf2{Intg_DR}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, ρ_ex = 0, growth = 1.15)
	
	i.dt[] == 0. && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)

	i.dt[] = min(i.dt[], sqrt(3) / ρ_ex)

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 DR stage 1:")
	if i.imex.im
		assign!(i.rhs, i.stage_1_rhs)
		assign!(i.w_2, i.rhs)
		
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.lhs((v, e, i.y[3]),)[1:2], i.w_2[2])
		newtonit!(elastic, i.w_2[1], i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage  = α -> i.rhs[3] - i.lhs((i.y[1], i.y[2], α),)[3]
		newtonit!(damage, i.w_2[3], i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_2, i.stage_1_rhs)
	end

	quiet || println("TR-BDF2 DR stage 2:")
    if i.imex.im
		assign!(i.rhs, i.stage_2_rhs)
		assign!(i.w_3, i.rhs)
		
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.lhs((v, e, i.y[3]),)[1:2], i.w_3[2])
		newtonit!(elastic, i.w_3[1], i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage  = α -> i.rhs[3] - i.lhs((i.y[1], i.y[2], α),)[3]
		newtonit!(damage,  i.w_3[3], i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_3, i.stage_2_rhs)
	end
	
	quiet || println("TR-BDF2 DR explicit error:")
	assign!(i.e_ex, expl_err)
	
	quiet || println("TR-BDF2 DR implicit error:")
	if i.imex.im
		assign!(i.rhs,  i.impl_err_rhs)
    	assign!(i.e_im, i.rhs)
    	
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.lhs((v, e, i.y[3]),)[1:2], i.e_im[2])
		newtonit!(elastic, i.e_im[1], i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage  = α -> i.rhs[3] - i.lhs((i.y[1], i.y[2], α),)[3]
		newtonit!(damage,  i.e_im[3], i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.e_im, i.impl_err_rhs)
	end

	quiet || println("TR-BDF2 DR update:")
	assign!(i.y, i.update)
	
	η_ex = map((e, y) -> l2(e) / (rtol * l2(y) + atol), i.e_ex, i.y)
	η_im = map((e, y) -> l2(e) / (rtol * l2(y) + atol), i.e_im, i.y)
	
	i.t[] += i.dt[]

    quiet || println("t = $(i.t[]), dt = $(i.dt[]), η_ex = $η_ex, η_im = $η_im")
	
	i.dt[] = min(1/sqrt(max(η_ex...)), 1/cbrt(max(η_im...)), growth) * i.dt[]

	return (η_ex, η_im)
end


function step!(i::tr_bdf2{Intg_DBR}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, ρ_ex = 0, growth = 1.15)
	
	i.dt[] == 0. && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)

	i.dt[] = min(i.dt[], sqrt(3) / ρ_ex)

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 DBR stage 1:")
	if i.imex.im
		assign!(i.rhs, i.stage_1_rhs)
		assign!(i.w_2, i.rhs)
		
		quiet || println("elasticity:")
		elastic  = SchurComplement(((v, e),) -> i.rhs[1:2] - i.lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_2.e)
		newtonit!(elastic,  i.w_2.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage   = α -> i.rhs[3] - i.lhs((i.y.v, i.y.e, α, i.y.β),)[3]
		newtonit!(damage,   i.w_2.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

		quiet || println("breakage:")
		breakage = β -> i.rhs[4] - i.lhs((i.y.v, i.y.e, i.y.α, β),)[4]
		newtonit!(breakage, i.w_2.β, i.h.β[1], i.h.β[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_2, i.stage_1_rhs)
	end

	quiet || println("TR-BDF2 DBR stage 2:")
    if i.imex.im
		assign!(i.rhs, i.stage_2_rhs)
		assign!(i.w_3, i.rhs)
		
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_3.e)
		newtonit!(elastic,  i.w_3.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage   = α -> i.rhs[3] - i.lhs((i.y.v, i.y.e, α, i.y.β),)[3]
		newtonit!(damage,   i.w_3.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("breakage:")
		breakage = β -> i.rhs[4] - i.lhs((i.y.v, i.y.e, i.y.α, β),)[4]
		newtonit!(breakage, i.w_3.β, i.h.β[1], i.h.β[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.w_3, i.stage_2_rhs)
	end
	
	quiet || println("TR-BDF2 DBR explicit error:")
	assign!(i.e_ex, expl_err)
	
	quiet || println("TR-BDF2 DBR implicit error:")
	if i.imex.im
		assign!(i.rhs,  i.impl_err_rhs)
    	assign!(i.e_im, i.rhs)
    	
		quiet || println("elasticity:")
		elastic  = SchurComplement(((v, e),) -> i.rhs[1:2] - i.lhs((v, e, i.y.α, i.y.β),)[1:2], i.e_im.e)
		newtonit!(elastic,  i.e_im.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage   = α -> i.rhs[3] - i.lhs((i.y.v, i.y.e, α, i.y.β),)[3]
		newtonit!(damage,   i.e_im.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("breakage:")
		breakage = β -> i.rhs[4] - i.lhs((i.y.v, i.y.e, i.y.α, β),)[4]
		newtonit!(breakage, i.e_im.β, i.h.β[1], i.h.β[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	else
		assign!(i.e_im, i.impl_err_rhs)
	end

	quiet || println("TR-BDF2 DBR update:")
	assign!(i.y, i.update)
	
	η_ex = map((e, y) -> l2(e) / (rtol * l2(y) + atol), i.e_ex, i.y)
	η_im = map((e, y) -> l2(e) / (rtol * l2(y) + atol), i.e_im, i.y)
	
	i.t[] += i.dt[]

    quiet || println("t = $(i.t[]), dt = $(i.dt[]), η_ex = $η_ex, η_im = $η_im")
	
	i.dt[] = min(1/sqrt(max(η_ex...)), 1/cbrt(max(η_im...)), growth) * i.dt[]

	return (η_ex, η_im)
end

