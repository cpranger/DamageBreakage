
abstract type Integrator end

struct tr_bdf2 <: Integrator
	f_ex
	f_im
    y
	w
	h
	r
	dw
	e
    h_t
end

tr_bdf2(f,          y; h_t = 1.) = tr_bdf2(x -> x, f,    y;  h_t = h_t)
tr_bdf2(f_ex, f_im, y; h_t     ) = tr_bdf2(f_ex,   f_im, y, [deepcopy(y) for _ in 1:3], [deepcopy(y) for _ in 1:5], deepcopy(y), deepcopy(y), deepcopy(y), h_t)

function step!(i::tr_bdf2; newton_maxit, newton_atol)
	assign!(i.w[1], i.y)
    
	stage_1 = w_2 -> (i.y + i.h_t * (
	    (2/1-sqrt(2)/1) * i.f_ex(i.w[1])
	  + (1/1-sqrt(2)/2) * i.f_im(i.w[1])
	)) - (w_2 - i.h_t * (
	    (1/1-sqrt(2)/2) * i.f_im(w_2)
	))

    println("TR-BDF2 stage 1:")
	newtonit!(stage_1, i.w[2], i.dw, i.r, i.h; maxit = newton_maxit, atol = newton_atol)

    stage_2 = w_3 -> (i.y + i.h_t * (
	    (1/2-sqrt(2)/3) * i.f_ex(i.w[1]) + (1/2+sqrt(2)/3) * i.f_ex(i.w[2])
	  + (    sqrt(2)/4) * i.f_im(i.w[1]) + (    sqrt(2)/4) * i.f_im(i.w[2])
	)) - (w_3 - i.h_t * (
	    (1/1-sqrt(2)/2) * i.f_im(w_3)
	))

    println("TR-BDF2 stage 2:")
    newtonit!(stage_2, i.w[3], i.dw, i.r, i.h; maxit = newton_maxit, atol = newton_atol)
	
	stage_e = e -> (i.h_t * (
	    ( 1/3-sqrt(2)/3) * i.f_im(i.w[1])
	  + ( 1/3          ) * i.f_im(i.w[2])
	  + (-2/3+sqrt(2)/3) * i.f_im(i.w[3])
	)) - (e - i.h_t * (
	    (1/1-sqrt(2)/2) * i.f_im(e)
	))

	println("TR-BDF2 error stage:")
    newtonit!(stage_e, i.e, i.dw, i.r, i.h; maxit = newton_maxit, atol = newton_atol)

    update = y -> y - (i.w[1] + i.h_t * (
        (    sqrt(2)/4)*(i.f_im(i.w[1]) + i.f_ex(i.w[1]))
      + (    sqrt(2)/4)*(i.f_im(i.w[2]) + i.f_ex(i.w[2]))
      + (1/1-sqrt(2)/2)*(i.f_im(i.w[3]) + i.f_ex(i.w[3]))
    ))
    
    println("TR-BDF2 update stage:")
    newtonit!(update, i.y, i.dw, i.r, i.h; maxit = newton_maxit, atol = newton_atol)

    return sqrt(dot(i.e, i.e))
end