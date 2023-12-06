
struct tr_bdf2
	f_ex
	f_im
	y
	w_1
	w_2
	w_3
	e
	h
	h_t
end

struct tr_bdf2_schur
	f_ex
	f_im
	y
	w_1
	w_2
	w_3
	e
	h
	h_t
end

tr_bdf2(      f,          y; h_t) = tr_bdf2(      y -> 0*y, f,    y;  h_t = h_t)
tr_bdf2(      f_ex, f_im, y; h_t) = tr_bdf2(      f_ex,     f_im, y, [deepcopy(y) for _ in 1:4]..., [deepcopy(y)    for _ in 1:7], h_t)

tr_bdf2_schur(f,          y; h_t) = tr_bdf2_schur(y -> 0*y, f,    y;  h_t = h_t)
tr_bdf2_schur(f_ex, f_im, y; h_t) = tr_bdf2_schur(f_ex,     f_im, y, [deepcopy(y) for _ in 1:4]..., [deepcopy(y[1]) for _ in 1:7], h_t)

tr_bdf2_stage_1(f_ex, f_im, w_1, h_t) = w_2 -> (w_1 + h_t * (
	(2/1-sqrt(2)/1) * f_ex(w_1)
  + (1/1-sqrt(2)/2) * f_im(w_1)
)) - (w_2 - h_t * (
	(1/1-sqrt(2)/2) * f_im(w_2)
))

tr_bdf2_stage_2(f_ex, f_im, w_1, w_2, h_t) = w_3 -> (w_1 + h_t * (
	(1/2-sqrt(2)/3) * f_ex(w_1) + (1/2+sqrt(2)/3) * f_ex(w_2)
  + (    sqrt(2)/4) * f_im(w_1) + (    sqrt(2)/4) * f_im(w_2)
)) - (w_3 - h_t * (
	(1/1-sqrt(2)/2) * f_im(w_3)
))

tr_bdf2_stage_e(f_ex, f_im, w_1, w_2, w_3, h_t) = e -> (h_t * (
	( 1/3-sqrt(2)/3) * f_im(w_1)
  + ( 1/3          ) * f_im(w_2)
  + (-2/3+sqrt(2)/3) * f_im(w_3)
)) - (e - h_t * (
	( 1/1-sqrt(2)/2) * f_im(e)
))

tr_bdf2_update(f_ex, f_im, w_1, w_2, w_3, h_t) = w_1 + h_t * (
	(    sqrt(2)/4)*(f_im(w_1) + f_ex(w_1))
  + (    sqrt(2)/4)*(f_im(w_2) + f_ex(w_2))
  + (1/1-sqrt(2)/2)*(f_im(w_3) + f_ex(w_3))
)

function init!(i::tr_bdf2, func; newton_maxit, newton_rtol)
	null = 0 * i.y
	assign!(i.y, func(null))
	newtonit!(func, i.y, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)
end

function init!(i::tr_bdf2_schur, func; newton_maxit, newton_rtol)
	null = 0 .* i.y
	assign!(i.y, func(null))
	sc = SchurComplement(func, i.y[2])
	newtonit!(sc, i.y[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)
end

function step!(i::tr_bdf2_schur; newton_maxit, newton_rtol)
	stage_1 = tr_bdf2_stage_1(i.f_ex, i.f_im, i.w_1, i.h_t)
	stage_2 = tr_bdf2_stage_2(i.f_ex, i.f_im, i.w_1, i.w_2, i.h_t)
	stage_e = tr_bdf2_stage_e(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.h_t)
	update  = tr_bdf2_update( i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.h_t)
	
	null = 0 * i.y

	assign!(i.w_1, i.y)

	println("TR-BDF2 stage 1:")
	assign!(i.w_2, stage_1(null))
	sc = SchurComplement(stage_1, i.w_2[2])
	newtonit!(sc, i.w_2[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)

    println("TR-BDF2 stage 2:")
    assign!(i.w_3, stage_2(null))
	sc = SchurComplement(stage_2, i.w_3[2])
	newtonit!(sc, i.w_3[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)

	println("TR-BDF2 error stage:")
	assign!(i.e,   stage_e(null))
    sc = SchurComplement(stage_e, i.e[2])
	newtonit!(sc, i.e[1],   i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)

	println("TR-BDF2 update:")
	assign!(i.y, update)

    return (l2(i.e[1]), l2(i.e[2]))
end

function step!(i::tr_bdf2; newton_maxit, newton_rtol)
	stage_1 = tr_bdf2_stage_1(i.f_ex, i.f_im, i.w_1, i.h_t)
	stage_2 = tr_bdf2_stage_2(i.f_ex, i.f_im, i.w_1, i.w_2, i.h_t)
	stage_e = tr_bdf2_stage_e(i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.h_t)
	update  = tr_bdf2_update( i.f_ex, i.f_im, i.w_1, i.w_2, i.w_3, i.h_t)
	
	null = 0 * i.y

	assign!(i.w_1, i.y)

	println("TR-BDF2 Schur stage 1:")
	assign!(i.w_2, stage_1(null))
	newtonit!(stage_1, i.w_2, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)

    println("TR-BDF2 Schur stage 2:")
    assign!(i.w_3, stage_2(null))
	newtonit!(stage_2, i.w_3, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)

	println("TR-BDF2 Schur error stage:")
	assign!(i.e,   stage_e(null))
    newtonit!(stage_e, i.e,   i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol)

	println("TR-BDF2 Schur update:")
	assign!(i.y, update)

    return l2(i.e)
end