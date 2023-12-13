
struct tr_bdf2
	f_ex
	f_im
	y
	w_1
	w_2
	w_3
	e
	h
	dt::Ref{Float64}
	 t::Ref{Float64}
end

struct tr_bdf2_schur
	f_ex
	f_im
	y
	w_1
	w_2
	w_3
	e
	h # copies of  u
	g # copies of (u, v)
	dt::Ref{Float64}
	 t::Ref{Float64}
end

tr_bdf2(      f,          y; dt = 0.) = tr_bdf2(      y -> 0*y, f,    y;  dt = dt)
tr_bdf2(      f_ex, f_im, y; dt = 0.) = tr_bdf2(      f_ex,     f_im, y, [deepcopy(y) for _ in 1:4]..., [deepcopy(y)    for _ in 1:7], Ref(dt), Ref(0.))

tr_bdf2_schur(f,          y; dt = 0.) = tr_bdf2_schur(y -> 0*y, f,    y;  dt = dt)
tr_bdf2_schur(f_ex, f_im, y; dt = 0.) = tr_bdf2_schur(f_ex,     f_im, y, [deepcopy(y) for _ in 1:4]..., [deepcopy(y[1]) for _ in 1:7], [deepcopy(y) for _ in 1:2], Ref(dt), Ref(0.))

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

tr_bdf2_stage_e(f_ex, f_im, w_1, w_2, w_3, dt) = e -> (dt * (
	( 1/3-sqrt(2)/3) * f_im(w_1)
  + ( 1/3          ) * f_im(w_2)
  + (-2/3+sqrt(2)/3) * f_im(w_3)
)) - (e - dt * (
	( 1/1-sqrt(2)/2) * f_im(e)
))

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

function step!(i::tr_bdf2_schur; newton_maxit, newton_rtol, quiet = false)
	stage_1 = tr_bdf2_stage_1(i.f_ex, i.f_im, i.w_1, i.dt[])
	stage_2 = tr_bdf2_stage_2(i.f_ex, i.f_im, i.w_1, i.w_2, i.dt[])
	stage_e = tr_bdf2_stage_e(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	update  = tr_bdf2_update( i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	
	null = 0 * i.y

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 Schur stage 1:")
	assign!(i.w_2, stage_1(null))
	sc = SchurComplement(stage_1, i.w_2[2])
	newtonit!(sc, i.w_2[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

    quiet || println("TR-BDF2 Schur stage 2:")
    assign!(i.w_3, stage_2(null))
	sc = SchurComplement(stage_2, i.w_3[2])
	newtonit!(sc, i.w_3[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 Schur error stage:")
	assign!(i.e,   stage_e(null))
    sc = SchurComplement(stage_e, i.e[2])
	newtonit!(sc, i.e[1],   i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 Schur update:")
	assign!(i.y, update)
	i.t[] += i.dt[]

    return (l2(i.e[1]), l2(i.e[2]))
end

function step!(i::tr_bdf2; newton_maxit, newton_rtol, quiet = false)
	stage_1 = tr_bdf2_stage_1(i.f_ex, i.f_im, i.w_1, i.dt[])
	stage_2 = tr_bdf2_stage_2(i.f_ex, i.f_im, i.w_1, i.w_2, i.dt[])
	stage_e = tr_bdf2_stage_e(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	update  = tr_bdf2_update( i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.dt[])
	
	null = 0 * i.y

	assign!(i.w_1, i.y)

	quiet || println("TR-BDF2 stage 1:")
	assign!(i.w_2, stage_1(null))
	newtonit!(stage_1, i.w_2, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 stage 2:")
    assign!(i.w_3, stage_2(null))
	newtonit!(stage_2, i.w_3, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 error stage:")
	assign!(i.e,   stage_e(null))
    newtonit!(stage_e, i.e,   i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

	quiet || println("TR-BDF2 update:")
	assign!(i.y, update)
	i.t[] += i.dt[]

    return l2(i.e)
end