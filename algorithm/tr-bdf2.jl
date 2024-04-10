
const nullfunc = y -> 0*y

abstract type IntgType end

struct Intg_    <: IntgType; end  # default for parabolic problems
struct Intg_SCR <: IntgType; end  # Schur Complement Reduction for hyperbolic problems
struct Intg_DR  <: IntgType; end  # Damage Rheology
struct Intg_DBR <: IntgType; end  # Damage-Breakage Rheology

struct tr_bdf2{T <: IntgType}
	f_im
	f_ex
	y
	e
	w_1
	w_2
	w_3
	w_4
	rhs
	h   #      helper vecs
	g   # more helper vecs
	tr1_lhs
	tr1_rhs
	bdf2_lhs
	bdf2_rhs
	tr2_lhs
	tr2_rhs
	update_y
	update_e
	imex
	t::FieldVal
	dt::FieldVal
end

tr_bdf2_data(                  y                              ) = (; y = y, e = deepcopy(y), w_1 = deepcopy(y), w_2 = deepcopy(y), w_3 = deepcopy(y), w_4 = deepcopy(y), rhs = deepcopy(y))
tr_bdf2_data(::Type{Intg_},    y                              ) = (; tr_bdf2_data(y)..., h = [deepcopy(y)    for _ in 1:7],                                                                            g = [])
tr_bdf2_data(::Type{Intg_SCR}, y::Tuple                       ) = (; tr_bdf2_data(y)..., h = [deepcopy(y[1]) for _ in 1:7],                                                                            g = [deepcopy(y)           for _ in 1:2])
tr_bdf2_data(::Type{Intg_DR},  y::NamedTuple{(:v, :e, :α)}    ) = (; tr_bdf2_data(values(y))..., h = (; v = [deepcopy(y.v) for _ in 1:7], α = [deepcopy(y.α) for _ in 1:7]),                                   g = [deepcopy((; e = y.e)) for _ in 1:2])
tr_bdf2_data(::Type{Intg_DBR}, y::NamedTuple{(:v, :e, :α, :β)}) = (; tr_bdf2_data(values(y))..., h = (; v = [deepcopy(y.v) for _ in 1:7], α = [deepcopy(y.α) for _ in 1:7], β = [deepcopy(y.β) for _ in 1:7]), g = [deepcopy((; e = y.e)) for _ in 1:2])

tr_bdf2_funcs(f_im, f_ex, data, dt::FieldVal) = (
	                                    tr1_lhs(      f_im,                                                 dt),
	prepare_assignment(data.w_1,        tr1_rhs(f_ex, f_im, data.w_1,                                       dt)),
	                                   bdf2_lhs(      f_im,                                                 dt),
	prepare_assignment(data.w_1,       bdf2_rhs(f_ex, f_im, data.w_1, data.w_2,                             dt)),
	                                    tr2_lhs(      f_im,                                                 dt),
	prepare_assignment(data.w_1,        tr2_rhs(f_ex, f_im, data.w_1, data.w_2,                             dt)),
	prepare_assignment(data.w_1, tr_bdf2_update(f_ex, f_im, data.w_1, data.w_2, data.w_3,                   dt)),
	prepare_assignment(data.w_1, tr_bdf2_error( f_ex, f_im, data.w_1,                     data.w_4, data.y, dt))
)

imex(f_im, f_ex) = (;
	im = typeof(f_im) != typeof(nullfunc),
	ex = typeof(f_ex) != typeof(nullfunc)
)

tr_bdf2(type, y; f_ex = nullfunc, f_im = nullfunc) = tr_bdf2_(;
	type = type,
	f_ex = f_ex,
	f_im = f_im, 
	data = tr_bdf2_data(type, y),
	dt   = Ref(Inf),
	t    = Ref(0.)
)
tr_bdf2_(; type, f_ex, f_im, data, t, dt) = tr_bdf2__(; 
	type  = type,
	f_ex  = f_ex,
	f_im  = f_im, 
	data  = data,
	funcs = tr_bdf2_funcs(f_im, f_ex, data, FieldVal(dt)),
	imex  = imex(f_im, f_ex),
	dt    = dt,
	t     = t
)
tr_bdf2__(      ; type, f_im, f_ex, data, funcs, imex, t, dt) = tr_bdf2{type}(f_im, f_ex, data..., funcs..., imex, FieldVal(t), FieldVal(dt))

global const sqrt2 = sqrt(2)

tr1_lhs(f_im, dt) = w_2 -> w_2 - dt * (1 - sqrt2/2) * f_im(w_2)

tr1_rhs(f_ex, f_im, w_1, dt) = w_1 + dt * (
	(2 - sqrt2  ) * f_ex(w_1)
  + (1 - sqrt2/2) * f_im(w_1)
)

bdf2_lhs(f_im, dt) = w_3 -> w_3 - dt * (1 - sqrt2/2) * f_im(w_3)

bdf2_rhs(f_ex, f_im, w_1, w_2, dt) = w_1 + dt * (
    (.5 - sqrt2/6) * f_ex(w_1) + sqrt2/4 * f_im(w_1)
  + (.5 + sqrt2/6) * f_ex(w_2) + sqrt2/4 * f_im(w_2)
)

tr2_lhs(f_im, dt) = w_4 -> w_4 - dt * .5 * f_im(w_4)

tr2_rhs(f_ex, f_im, w_1, w_2, dt) = w_1 + dt * (
    (5/6 - sqrt2/12) * f_ex(w_1) + .5 * f_im(w_1)
  + (1/6 + sqrt2/12) * f_ex(w_2)
)

tr_bdf2_update(f_ex, f_im, w_1, w_2, w_3, dt) = w_1 + dt * (
    (    sqrt2/4)*(f_im(w_1) + f_ex(w_1))
  + (    sqrt2/4)*(f_im(w_2) + f_ex(w_2))
  + (1 - sqrt2/2)*(f_im(w_3) + f_ex(w_3))
)

# to be computed after update
tr_bdf2_error(f_ex, f_im, w_1, w_4, y, dt) = (12 + 8sqrt2) * (y - w_1 - dt * (
    .5*(f_im(w_1) + f_ex(w_1))
  + .5*(f_im(w_4) + f_ex(w_4))
))


# err_lhs(f_im, dt) = w -> w - dt * .5 * f_im(w)

# err_rhs(f_ex, f_im, w_1, w_2, w_3, dt) = w_1 + dt * (
# 	(5/6 - sqrt2/12) * f_ex(w_1) + .5 * f_im(w_1) +
# 	(1/6 + sqrt2/12) * f_ex(w_2) - w_3
# )


# tr_bdf2_lhs(f_im, dt) = w -> w - dt * (1/1-sqrt(2)/2) * f_im(w)

# tr_bdf2_stage_2_rhs(f_ex, f_im, y, w_1, w_2, dt) = y + dt * (
# 	(1/2-sqrt(2)/3) * f_ex(w_1) + (sqrt(2)/4) * f_im(w_1)
#   + (1/2+sqrt(2)/3) * f_ex(w_2) + (sqrt(2)/4) * f_im(w_2)
# )

# tr_bdf2_impl_err_rhs(f_im, w_1, w_2, w_3, dt) = dt * (
# 	( 1/3-sqrt(2)/3) * f_im(w_1)
#   + ( 1/3          ) * f_im(w_2)
#   + (-2/3+sqrt(2)/3) * f_im(w_3)
# )

# tr_bdf2_expl_err(f_ex, w_1, w_2, w_3, dt) = dt * (
# 	( 1/2-sqrt(2)/2) * f_ex(w_1)
#   + ( 1/2          ) * f_ex(w_2)
#   + (-1/1+sqrt(2)/2) * f_ex(w_3)
# )

# tr_bdf2_update(f_ex, f_im, w_1, w_2, w_3, dt) = w_1 + dt * (
# 	(    sqrt(2)/4)*(f_im(w_1) + f_ex(w_1))
#   + (    sqrt(2)/4)*(f_im(w_2) + f_ex(w_2))
#   + (1/1-sqrt(2)/2)*(f_im(w_3) + f_ex(w_3))
# )

function init!(i::tr_bdf2{Intg_}, func; newton_maxit, newton_rtol, quiet = false)
	null = 0 * i.y
	assign!(i.y, func(null))
	newtonit!(func, i.y, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
end

function init_stepsize!(i::tr_bdf2; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, ρ_ex = Inf)
	dt = i.dt[]
	t  = i.t[]
	for _ in 1:10
		i.dt[] = dt
		η = step!(i,
			rtol = rtol,
			atol = atol,
			newton_maxit = newton_maxit,
			newton_rtol = newton_rtol,
			quiet = quiet,
			ρ_ex = ρ_ex
		)
		i.t[]  = t
		assign!(i.y, i.w_1)
		#=quiet ||=# println("dt = $dt: η = $η")
		abs(max(η...) - 1) < .1 && break
		dt /= sqrt(max(η...))
	end

	i.dt[] = dt
	return dt
end


function step!(i::tr_bdf2{Intg_}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, growth = 1.15, ρ_ex = Inf)
	
	do_init = (i.dt[] == Inf)
	
	# TODO: improve stability criterion
	if ρ_ex == Inf
		circles = gershgorin!(linearize(i.f_ex, i.y), (i.w_1, i.w_2, i.w_3))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end

	i.dt[] = min(i.dt[], 3/sqrt2/ρ_ex)
	
	do_init && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)


	assign!(i.w_1, i.y)

	quiet || println("first trapezoidal stage:")
	assign!(i.w_2, i.tr1_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_2)
		newtonit!(w -> i.rhs - i.tr1_lhs(w), i.w_2, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("bdf2 stage:")
	assign!(i.w_3, i.bdf2_rhs)
    if i.imex.im
		assign!(i.rhs, i.w_3)
		newtonit!(w -> i.rhs - i.bdf2_lhs(w), i.w_3, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("second trapezoidal stage:")
	assign!(i.w_4, i.tr2_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_4)
		newtonit!(w -> i.rhs - i.tr2_lhs(w), i.w_4, i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("finishing stage:")
	assign!(i.y, i.update_y)
	
	quiet || println("error stage:")
	assign!(i.e, i.update_e)
    
	i.t[] += i.dt[]

	η = l2(i.e) / (rtol * l2(i.y) + atol)

	quiet || println("t = $(i.t[]), dt = $(i.dt[]), η = $η")
	
	i.dt[] = min(1/sqrt(max(η...)), growth) * i.dt[]

	return η
end


function init!(i::tr_bdf2{Intg_SCR}, func; newton_maxit, newton_rtol, quiet = false)
	null = 0 .* i.y
	assign!(i.y, func(null))
	sc = SchurComplement(func, i.y[2])
	newtonit!(sc, i.y[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
end

function step!(i::tr_bdf2{Intg_SCR}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, growth = 1.15, ρ_ex = Inf)
	
	do_init = (i.dt[] == Inf)
	
	# TODO: improve stability criterion
	if ρ_ex == Inf
		circles = gershgorin!(linearize(i.f_ex, i.y), (i.w_1, i.w_2, i.w_3))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end

	i.dt[] = min(i.dt[], 3/sqrt2/ρ_ex)
	
	do_init && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)

	assign!(i.w_1, i.y)

	quiet || println("first trapezoidal stage:")
	assign!(i.w_2, i.tr1_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_2)
		sc = SchurComplement(w -> i.rhs - i.tr1_lhs(w), i.w_2[2])
		newtonit!(sc, i.w_2[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

    quiet || println("bdf2 stage:")
	assign!(i.w_3, i.bdf2_rhs)
    if i.imex.im
		assign!(i.rhs, i.w_3)
		sc = SchurComplement(w -> i.rhs - i.bdf2_lhs(w), i.w_3[2])
		newtonit!(sc, i.w_3[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("second trapezoidal stage:")
	assign!(i.w_4, i.tr2_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_4)
		sc = SchurComplement(w -> i.rhs - i.tr2_lhs(w), i.w_4[2])
		newtonit!(sc, i.w_4[1], i.h[1], i.h[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("finishing stage:")
	assign!(i.y, i.update_y)
	
	quiet || println("error stage:")
	assign!(i.e, i.update_e)
    
	i.t[] += i.dt[]

	η = l2(i.e) / (rtol * l2(i.y) + atol)

	quiet || println("t = $(i.t[]), dt = $(i.dt[]), η = $η")
	
	i.dt[] = min(1/sqrt(max(η...)), growth) * i.dt[]

	return η
end

function step!(i::tr_bdf2{Intg_DR}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, growth = 1.15, ρ_ex = Inf)
	
	do_init = (i.dt[] == Inf)
	
	# TODO: improve stability criterion
	if ρ_ex == Inf
		circles = gershgorin!(linearize(((e, α),) -> i.f_ex((i.y[1], e, α),)[2:3], i.y[2:3]), (i.w_1[2:3], i.w_2[2:3], i.w_3[2:3]))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end

	i.dt[] = min(i.dt[], 3/sqrt2/ρ_ex)
	
	do_init && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)

	assign!(i.w_1, i.y)

	quiet || println("first trapezoidal stage:")
	assign!(i.w_2, i.tr1_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_2)
		
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.tr1_lhs((v, e, i.y[3]),)[1:2], i.w_2[2])
		newtonit!(elastic, i.w_2[1], i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage  = α -> i.rhs[3] - i.tr1_lhs((i.y[1], i.y[2], α),)[3]
		newtonit!(damage, i.w_2[3], i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("bdf2 stage:")
    assign!(i.w_3, i.bdf2_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_3)
		
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.bdf2_lhs((v, e, i.y[3]),)[1:2], i.w_3[2])
		newtonit!(elastic, i.w_3[1], i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage  = α -> i.rhs[3] - i.bdf2_lhs((i.y[1], i.y[2], α),)[3]
		newtonit!(damage,  i.w_3[3], i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end
	
	quiet || println("second trapezoidal stage:")
	assign!(i.w_4, i.tr2_rhs)
	if i.imex.im
		assign!(i.rhs,  i.w_4)
    	
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.tr2_lhs((v, e, i.y[3]),)[1:2], i.w_4[2])
		newtonit!(elastic, i.w_4[1], i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage  = α -> i.rhs[3] - i.tr2_lhs((i.y[1], i.y[2], α),)[3]
		newtonit!(damage,  i.w_4[3], i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("finishing stage:")
	assign!(i.y, i.update_y)
	
	quiet || println("error stage:")
	assign!(i.e, i.update_e)
    
	i.t[] += i.dt[]

	η = l2(i.e) / (rtol * l2(i.y) + atol)

	quiet || println("t = $(i.t[]), dt = $(i.dt[]), η = $η")
	
	i.dt[] = min(1/sqrt(max(η...)), growth) * i.dt[]

	return η
end


function step!(i::tr_bdf2{Intg_DBR}; rtol, atol = 0, newton_maxit, newton_rtol, quiet = false, growth = 1.15, ρ_ex = Inf)
	
	do_init = (i.dt[] == Inf)
	
	# TODO: improve stability criterion
	if ρ_ex == Inf
		circles = gershgorin!(linearize(((e, α, β),) -> i.f_ex((i.y[1], e, α, β),)[2:4], i.y[2:4]), (i.w_1[2:4], i.w_2[2:4], i.w_3[2:4]))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end

	i.dt[] = min(i.dt[], 3/sqrt2/ρ_ex)
	
	do_init && init_stepsize!(i;
		rtol = rtol,
		atol = atol,
		newton_maxit = newton_maxit,
		newton_rtol = newton_rtol,
		quiet = quiet,
		ρ_ex = ρ_ex
	)

	assign!(i.w_1, i.y)

	quiet || println("first trapezoidal stage:")
	assign!(i.w_2, i.tr1_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_2)
		
		quiet || println("elasticity:")
		elastic  = SchurComplement(((v, e),) -> i.rhs[1:2] - i.tr1_lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_2.e)
		newtonit!(elastic,  i.w_2.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage   = α -> i.rhs[3] - i.tr1_lhs((i.y.v, i.y.e, α, i.y.β),)[3]
		newtonit!(damage,   i.w_2.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)

		quiet || println("breakage:")
		breakage = β -> i.rhs[4] - i.tr1_lhs((i.y.v, i.y.e, i.y.α, β),)[4]
		newtonit!(breakage, i.w_2.β, i.h.β[1], i.h.β[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("bdf2 stage:")
    assign!(i.w_3, i.bdf2_rhs)
	if i.imex.im
		assign!(i.rhs, i.w_3)
		
		quiet || println("elasticity:")
		elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.bdf2_lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_3.e)
		newtonit!(elastic,  i.w_3.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage   = α -> i.rhs[3] - i.bdf2_lhs((i.y.v, i.y.e, α, i.y.β),)[3]
		newtonit!(damage,   i.w_3.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("breakage:")
		breakage = β -> i.rhs[4] - i.bdf2_lhs((i.y.v, i.y.e, i.y.α, β),)[4]
		newtonit!(breakage, i.w_3.β, i.h.β[1], i.h.β[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end
	
	quiet || println("second trapezoidal stage:")
	assign!(i.w_4, i.tr2_rhs)
	if i.imex.im
		assign!(i.rhs,  i.w_4)
    	
		quiet || println("elasticity:")
		elastic  = SchurComplement(((v, e),) -> i.rhs[1:2] - i.tr2_lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_4.e)
		newtonit!(elastic,  i.w_4.v, i.h.v[1], i.h.v[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("damage:")
		damage   = α -> i.rhs[3] - i.tr2_lhs((i.y.v, i.y.e, α, i.y.β),)[3]
		newtonit!(damage,   i.w_4.α, i.h.α[1], i.h.α[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
		
		quiet || println("breakage:")
		breakage = β -> i.rhs[4] - i.tr2_lhs((i.y.v, i.y.e, i.y.α, β),)[4]
		newtonit!(breakage, i.w_4.β, i.h.β[1], i.h.β[2:end]; maxit = newton_maxit, rtol = newton_rtol, quiet = quiet)
	end

	quiet || println("finishing stage:")
	assign!(i.y, i.update_y)
	
	quiet || println("error stage:")
	assign!(i.e, i.update_e)
    
	i.t[] += i.dt[]

	η = l2(i.e) / (rtol * l2(i.y) + atol)

	quiet || println("t = $(i.t[]), dt = $(i.dt[]), η = $η")
	
	i.dt[] = min(1/sqrt(max(η...)), growth) * i.dt[]

	return η
end

