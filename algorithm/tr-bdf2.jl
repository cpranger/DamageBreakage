
struct tr_bdf2
	f_ex
	f_im
	y
	w_1
	w_2
	w_3
	e_im
	e_ex
	h
	t::Ref{Float64}
	dt::Ref{Float64}
end

struct tr_bdf2_schur
	f_ex
	f_im
	y
	w_1
	w_2
	w_3
	e_im
	e_ex
	h # copies of  u
	g # copies of (u, v)
	t::Ref{Float64}
	dt::Ref{Float64}
end

struct tr_bdf2_dr
	f_ex
	f_im
	y
	w_1
	w_2
	w_3
	e_im
	e_ex
	h # copies of (v,  , α) for implicit solves
	g # copies of (   e,  ) for other purposes
	t::Ref{Float64}
	dt::Ref{Float64}
end

tr_bdf2(      f,          y)               = tr_bdf2(      y -> 0*y, f,    y)
tr_bdf2(      f,          y; dt::Ref)      = tr_bdf2(      y -> 0*y, f,    y;  dt = dt)
tr_bdf2(      f_ex, f_im, y; dt = Ref(0.)) = tr_bdf2(      f_ex,     f_im, y, [deepcopy(y) for _ in 1:5]..., [deepcopy(y)    for _ in 1:7], Ref(0.), dt)

tr_bdf2_schur(f,          y)               = tr_bdf2_schur(y -> 0*y, f,    y)
tr_bdf2_schur(f,          y; dt::Ref)      = tr_bdf2_schur(y -> 0*y, f,    y;  dt = dt)
tr_bdf2_schur(f_ex, f_im, y; dt = Ref(0.)) = tr_bdf2_schur(f_ex,     f_im, y, [deepcopy(y) for _ in 1:5]..., [deepcopy(y[1]) for _ in 1:7], [deepcopy(y) for _ in 1:2], Ref(0.), dt)

tr_bdf2_dr(   f_ex, f_im, y::Tuple;      dt = Ref(0.)) = tr_bdf2_dr(f_ex, f_im, NamedTuple{(:v, :e, :α)}(y); dt = dt)
tr_bdf2_dr(   f_ex, f_im, y::NamedTuple; dt = Ref(0.)) = tr_bdf2_dr(f_ex, f_im, y, [deepcopy(y) for _ in 1:5]..., (; v = [deepcopy(y.v) for _ in 1:7], α = [deepcopy(y.α) for _ in 1:7]), [deepcopy((; e = y.e)) for _ in 1:2], Ref(0.), dt)

tr_bdf2_stage_1(f_ex, f_im, w_1, dt) = w_2 -> (w_1 + dt * (
	(2/1-sqrt(2)/1) * f_ex(w_1)
  + (1/1-sqrt(2)/2) * f_im(w_1)
)) - (w_2 - dt * (
	(1/1-sqrt(2)/2) * f_im(w_2)
))

tr_bdf2_stage_2(f_ex, f_im, w_1, w_2, dt) = w_3 -> (w_1 + dt * (
	(1/2-sqrt(2)/3) * f_ex(w_1) + (1/2+sqrt(2)/3) * f_ex(w_2)
  + (    sqrt(2)/4) * f_im(w_1) + (    sqrt(2)/4) * f_im(w_2)
)) - (w_3 - dt * (
	(1/1-sqrt(2)/2) * f_im(w_3)
))

tr_bdf2_impl_err(f_ex, f_im, w_1, w_2, w_3, dt) = e -> (dt * (
	( 1/3-sqrt(2)/3) * f_im(w_1)
  + ( 1/3          ) * f_im(w_2)
  + (-2/3+sqrt(2)/3) * f_im(w_3)
)) - (e - dt * (
	( 1/1-sqrt(2)/2) * f_im(e)
))

tr_bdf2_expl_err(f_ex, f_im, w_1, w_2, w_3, dt) = dt * (
	( 1/2-sqrt(2)/2) * f_ex(w_1)
  + ( 1/2          ) * f_ex(w_2)
  + (-1/1+sqrt(2)/2) * f_ex(w_3)
)

tr_bdf2_update(f_ex, f_im, w_1, w_2, w_3, dt) = w_1 + dt * (
	(    sqrt(2)/4)*(f_im(w_1) + f_ex(w_1))
  + (    sqrt(2)/4)*(f_im(w_2) + f_ex(w_2))
  + (1/1-sqrt(2)/2)*(f_im(w_3) + f_ex(w_3))
)

function init!(i::tr_bdf2, func; newton_maxit, newton_rtol, quiet = false)
	null = 0 * i.y
	assign!(i.y, func(null))
	newtonit!(func, i.y, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
end

function init!(i::tr_bdf2_schur, func; newton_maxit, newton_rtol, quiet = false)
	null = 0 .* i.y
	assign!(i.y, func(null))
	sc = SchurComplement(func, i.y[2])
	newtonit!(sc, i.y[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
end

function step!(i::tr_bdf2_schur; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false)
	stage_1  = tr_bdf2_stage_1( i.f_ex, i.f_im, i.w_1, i.dt[])
	stage_2  = tr_bdf2_stage_2( i.f_ex, i.f_im, i.w_1, i.w_2, i.dt[])
	impl_err = tr_bdf2_impl_err(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	expl_err = tr_bdf2_expl_err(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	update   = tr_bdf2_update(  i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	
	null = 0 * i.y

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 Schur stage 1:")
	assign!(i.w_2,   stage_1(null))
	sc = SchurComplement(stage_1, i.w_2[2])
	newtonit!(sc, i.w_2[1],  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

    quiet || println("TR-BDF2 Schur stage 2:")
    assign!(i.w_3,   stage_2(null))
	sc = SchurComplement(stage_2, i.w_3[2])
	newtonit!(sc, i.w_3[1],  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 Schur explicit error:")
	assign!(i.e_ex, expl_err)
	
	quiet || println("TR-BDF2 Schur implicit error:")
	assign!(i.e_im, impl_err(null))
    sc = SchurComplement(impl_err, i.e_im[2])
	newtonit!(sc, i.e_im[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 Schur update:")
	assign!(i.y, update)
	
	i.t[] += i.dt[]

    η_ex = l2(i.e_ex) / (rtol * l2(i.y) + atol)
	η_im = l2(i.e_im) / (rtol * l2(i.y) + atol)

	return (η_ex, η_im)
end





function step!(i::tr_bdf2_dr; rtol, atol = 0, newton_maxit, newton_atol = 0, newton_rtol, quiet = false, ρ_ex = 0, safety = .95, growth = 1.15)
	stage_1  = tr_bdf2_stage_1( i.f_ex, i.f_im, values(i.w_1), i.dt[])
	stage_2  = tr_bdf2_stage_2( i.f_ex, i.f_im, values(i.w_1), values(i.w_2), i.dt[])
	impl_err = tr_bdf2_impl_err(i.f_ex, i.f_im, values(i.w_1), values(i.w_2), values(i.w_3), i.dt[])
	expl_err = tr_bdf2_expl_err(i.f_ex, i.f_im, values(i.w_1), values(i.w_2), values(i.w_3), i.dt[])
	update   = tr_bdf2_update(  i.f_ex, i.f_im, values(i.w_1), values(i.w_2), values(i.w_3), i.dt[])
	
	i.dt[] = min(i.dt[], safety * sqrt(3) / ρ_ex)

	null = 0 * i.y

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 DR stage 1:")
	assign!(i.w_2,   stage_1(null))
	ve = SchurComplement(((v, e),) -> stage_1((v, e, null.α),)[1:2], i.w_2.e)
	α  = α -> stage_1((null.v, null.e, α),)[3]
	newtonit!(ve, i.w_2.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	newtonit!(α,  i.w_2.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

    quiet || println("TR-BDF2 DR stage 2:")
    assign!(i.w_3,   stage_2(null))
	ve = SchurComplement(((v, e),) -> stage_2((v, e, null.α),)[1:2], i.w_3.e)
	α  = α -> stage_2((null.v, null.e, α),)[3]
	newtonit!(ve, i.w_3.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	newtonit!(α,  i.w_3.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	
	quiet || println("TR-BDF2 DR explicit error:")
	assign!(i.e_ex, expl_err)
	
	quiet || println("TR-BDF2 DR implicit error:")
	assign!(i.e_im, impl_err(null))
    ve = SchurComplement(((v, e),) -> impl_err((v, e, null.α),)[1:2], i.e_im.e)
	α  = α -> impl_err((null.v, null.e, α),)[3]
	newtonit!(ve, i.e_im.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	newtonit!(α,  i.e_im.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 DR update:")
	assign!(i.y, update)
	
	η_ex = map((e, y) -> l2(e) / (rtol * l2(y) + atol), i.e_ex, i.y)
	η_ex = map((e, y) -> l2(e) / (rtol * l2(y) + atol), i.e_im, i.y)
	
	quiet || println("t = $(i.t[]), dt = $(i.dt[]), η_ex = $η_ex, η_im = $η_im")
	
	i.t[] += i.dt[]

    i.dt[] = min(1/sqrt(max(η_ex...)), 1/cbrt(max(η_im...)), growth) * i.dt[]

	return (i.t[], η_ex, η_im)
end







function step!(i::tr_bdf2; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false)
	stage_1  = tr_bdf2_stage_1( i.f_ex, i.f_im, i.w_1, i.dt[])
	stage_2  = tr_bdf2_stage_2( i.f_ex, i.f_im, i.w_1, i.w_2, i.dt[])
	impl_err = tr_bdf2_impl_err(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	expl_err = tr_bdf2_expl_err(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	update   = tr_bdf2_update(  i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	
	null = 0 * i.y

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 stage 1:")
	assign!(i.w_2, stage_1(null))
	newtonit!(stage_1,  i.w_2,  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 stage 2:")
    assign!(i.w_3, stage_2(null))
	newtonit!(stage_2,  i.w_3,  i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 explicit error:")
	assign!(i.e_ex, expl_err)
    
	quiet || println("TR-BDF2 implicit error:")
	assign!(i.e_im, impl_err(null))
    newtonit!(impl_err, i.e_im, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 update:")
	assign!(i.y, update)
	
	i.t[] += i.dt[]

	η_ex = max(map((u, ε) -> l2(ε) / (rtol * l2(u) + atol), i.y, i.e_ex)...)
	η_im = max(map((u, ε) -> l2(ε) / (rtol * l2(u) + atol), i.y, i.e_im)...)

	return (i.t[], η_ex, η_im)
end