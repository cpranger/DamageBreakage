
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
	err_im
	err_ex
	imex
	t::FieldVal
	dt::FieldVal
	dt_::FieldVal
	k::Ref{Int}
end

tr_bdf2_data(                  y                              ) = (; y = y, e = deepcopy(y), w_1 = deepcopy(y), w_2 = deepcopy(y), w_3 = deepcopy(y), w_4 = deepcopy(y), rhs = deepcopy(y))
tr_bdf2_data(::Type{Intg_},    y                              ) = (; tr_bdf2_data(y)...,         h = [deepcopy(y)    for _ in 1:7],                                                                            g = [])
tr_bdf2_data(::Type{Intg_SCR}, y::Tuple                       ) = (; tr_bdf2_data(y)...,         h = [deepcopy(y[1]) for _ in 1:7],                                                                            g = [deepcopy(y)           for _ in 1:2])
tr_bdf2_data(::Type{Intg_DR},  y::NamedTuple{(:u, :v, :e, :α)}) = (; tr_bdf2_data(values(y))..., h = (; v = [deepcopy(y.v) for _ in 1:7], α = [deepcopy(y.α) for _ in 1:7]),                                   g = [deepcopy((; e = y.e)) for _ in 1:2])
tr_bdf2_data(::Type{Intg_DBR}, y::NamedTuple{(:v, :e, :α, :β)}) = (; tr_bdf2_data(values(y))..., h = (; v = [deepcopy(y.v) for _ in 1:7], α = [deepcopy(y.α) for _ in 1:7], β = [deepcopy(y.β) for _ in 1:7]), g = [deepcopy((; e = y.e)) for _ in 1:2])

tr_bdf2_funcs(f_im, f_ex, data, t::FieldVal, dt::FieldVal) = (
	                                    tr1_lhs(      f_im,                                         t, dt),
	prepare_assignment(data.w_1,        tr1_rhs(f_ex, f_im, data.w_1,                               t, dt)),
	                                   bdf2_lhs(      f_im,                                         t, dt),
	prepare_assignment(data.w_1,       bdf2_rhs(f_ex, f_im, data.w_1, data.w_2,                     t, dt)),
	                                    tr2_lhs(      f_im,                                         t, dt),
	prepare_assignment(data.w_1,        tr2_rhs(f_ex, f_im, data.w_1, data.w_2,                     t, dt)),
	prepare_assignment(data.w_1, tr_bdf2_update(f_ex, f_im, data.w_1, data.w_2, data.w_3,           t, dt)),
	prepare_assignment(data.w_1,  tr_bdf2_error(      f_im, data.w_1, data.w_2, data.w_3, data.w_4, t, dt)),
	prepare_assignment(data.w_1,  tr_bdf2_error(f_ex,       data.w_1, data.w_2, data.w_3, data.w_4, t, dt))
)

imex(f_im, f_ex) = (;
	im = typeof(f_im) != typeof(nullfunc),
	ex = typeof(f_ex) != typeof(nullfunc)
)

tr_bdf2(type, y; f_ex = nullfunc, f_im = nullfunc, dt = Ref(Inf), dt_ = Ref(0.), t = Ref(BigFloat(0.)), k = Ref(0)) = tr_bdf2_(;
	type = type,
	f_ex = f_ex,
	f_im = f_im, 
	data = tr_bdf2_data(type, y),
	k    = k,
	dt   = dt,
	dt_  = dt_,
	t    = t
)
tr_bdf2_(; type, f_ex, f_im, data, t, dt, dt_, k) = tr_bdf2__(; 
	type  = type,
	f_ex  = f_ex,
	f_im  = f_im, 
	data  = data,
	funcs = tr_bdf2_funcs(f_im, f_ex, data, FieldVal(t), FieldVal(dt)),
	imex  = imex(f_im, f_ex),
	k     = k,
	dt    = dt,
	dt_   = dt_,
	t     = t
)
tr_bdf2__(      ; type, f_im, f_ex, data, funcs, imex, t, dt, dt_, k) = tr_bdf2{type}(f_im, f_ex, data..., funcs..., imex, FieldVal(t), FieldVal(dt), FieldVal(dt_), k)

global const sqrt2 = sqrt(2)

tr1_lhs(f_im, t, dt) = w_2 -> w_2 - dt * (1 - sqrt2/2) * f_im(t + (2 - sqrt2)*dt, w_2)

tr1_rhs(f_ex, f_im, w_1, t, dt) = w_1 + dt * (
	(2 - sqrt2  ) * f_ex(t, w_1)
  + (1 - sqrt2/2) * f_im(t, w_1)
)

bdf2_lhs(f_im, t, dt) = w_3 -> w_3 - dt * (1 - sqrt2/2) * f_im(t + dt, w_3)

bdf2_rhs(f_ex, f_im, w_1, w_2, t, dt) = w_1 + dt * (
    (.5 - sqrt2/6) * f_ex(t                 , w_1) + sqrt2/4 * f_im(t                 , w_1)
  + (.5 + sqrt2/6) * f_ex(t + (2 - sqrt2)*dt, w_2) + sqrt2/4 * f_im(t + (2 - sqrt2)*dt, w_2)
)

tr2_lhs(f_im, t, dt) = w_4 -> w_4 - dt * .5 * f_im(t + dt, w_4)

tr2_rhs(f_ex, f_im, w_1, w_2, t, dt) = w_1 + dt * (
    (5/6 - sqrt2/12) * f_ex(t                 , w_1) + .5 * f_im(t, w_1)
  + (1/6 + sqrt2/12) * f_ex(t + (2 - sqrt2)*dt, w_2)
)

tr_bdf2_update(f_ex, f_im, w_1, w_2, w_3, t, dt) = w_1 + dt * (
    (    sqrt2/4) * (f_im(t                 , w_1) + f_ex(t                 , w_1))
  + (    sqrt2/4) * (f_im(t + (2 - sqrt2)*dt, w_2) + f_ex(t + (2 - sqrt2)*dt, w_2))
  + (1 - sqrt2/2) * (f_im(t +             dt, w_3) + f_ex(t +             dt, w_3))
)

tr_bdf2_error(f, w_1, w_2, w_3, w_4, t, dt) = (12 + 8sqrt2) * dt * (
    (-.5 + sqrt2/4) * f(t                 , w_1)
  + (      sqrt2/4) * f(t + (2 - sqrt2)*dt, w_2)
  + (  1 - sqrt2/2) * f(t +             dt, w_3)
  + (-.5          ) * f(t +             dt, w_4)
)

function init!(i::tr_bdf2{Intg_}, func; newton_maxit, newton_rtol)
	null = 0 * i.y
	assign!(i.y, func(null))
	newtonit!(
		func, i.y, i.h[1], i.h[2:end];
		maxit    = newton_maxit,
		rtol     = newton_rtol,
		cg_maxit = cg_maxit,
		cg_rtol  = cg_rtol
	)
end

function step!(i::tr_bdf2{Intg_}; rtol, atol = 0, newton_maxit, newton_rtol, cg_maxit, cg_rtol, growth, safety, ρ_ex = Inf)
	
	if ρ_ex == 0 && i.dt[] == Inf # fully implicit
		i.dt[] = 1 # some value, to be calibrated
	end

	# TODO: improve stability criterion
	if ρ_ex == Inf
		circles = gershgorin!(linearize(y -> i.f_ex(i.t[], y), i.y), (i.w_1, i.w_2, i.w_3))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end

	i.dt_[] = safety * (3/sqrt2) / ρ_ex
	i.dt[]  = min(i.dt[], i.dt_[])
	
	assign!(i.w_1, i.y)

	@algo_step begin
		@verbo_println("first trapezoidal stage:")
		assign!(i.w_2, i.tr1_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_2)
			newtonit!(
				w -> i.rhs - i.tr1_lhs(w), i.w_2, i.h[1], i.h[2:end];
				maxit    = newton_maxit,
				rtol     = newton_rtol,
				cg_maxit = cg_maxit,
				cg_rtol  = cg_rtol
			)
		end
	end

	@algo_step begin
		@verbo_println("bdf2 stage:")
		assign!(i.w_3, i.bdf2_rhs)
    	if i.imex.im
			assign!(i.rhs, i.w_3)
			newtonit!(
				w -> i.rhs - i.bdf2_lhs(w), i.w_3, i.h[1], i.h[2:end];
				maxit    = newton_maxit,
				rtol     = newton_rtol,
				cg_maxit = cg_maxit,
				cg_rtol  = cg_rtol
			)
		end
	end

	@algo_step begin
		@verbo_println("second trapezoidal stage:")
		assign!(i.w_4, i.tr2_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_4)
			newtonit!(
				w -> i.rhs - i.tr2_lhs(w), i.w_4, i.h[1], i.h[2:end];
				maxit    = newton_maxit,
				rtol     = newton_rtol,
				cg_maxit = cg_maxit,
				cg_rtol  = cg_rtol
			)
		end
	end

	assign!(i.e, i.err_ex + i.err_im)
	
	# note: relative to w_3 instead of new y
	η = l2(i.e) / (rtol * l2(i.w_3) + atol)
	
	@verbo_println("TR-BDF2 $(i.k[]), t = $(@sprintf "%.3e" i.t[]), dt = $(@sprintf "%.3e" i.dt[]), η = $(@sprintf "%.3e" η)")
	
	if i.k[] == 0 && i.dt[] < i.dt_[] && abs(η - 1) > .1
		i.dt[] *= 0.1 + 0.9/sqrt(η)
		return step!(i;
			rtol         = rtol,
			atol         = atol,
			newton_maxit = newton_maxit,
			newton_rtol  = newton_rtol,
			cg_maxit     = cg_maxit,
			cg_rtol      = cg_rtol,
			growth       = growth,
			safety       = safety,
			ρ_ex         = ρ_ex
		)
	end

	assign!(i.y, i.update_y)
	
	i.t[] += i.dt[]
	i.k[] += 1

	i.dt[] = min(1/sqrt(η), growth) * i.dt[]

	return η
end

function init!(i::tr_bdf2{Intg_SCR}, func; newton_maxit, newton_rtol, cg_maxit, cg_rtol)
	null = 0 .* i.y
	assign!(i.y, func(null))
	sc = SchurComplement(func, i.y[2])
	newtonit!(
		sc, i.y[1], i.h[1], i.h[2:end];
		maxit    = newton_maxit,
		rtol     = newton_rtol,
		cg_maxit = cg_maxit,
		cg_rtol  = cg_rtol
	)
end

function step!(i::tr_bdf2{Intg_SCR}; rtol, atol = 0, newton_maxit, newton_rtol, cg_maxit, cg_rtol, growth, safety, ρ_ex = Inf)
	
	if ρ_ex == 0 && i.dt[] == Inf # fully implicit
		i.dt[] = 1 # some value, to be calibrated
	end

	# TODO: improve stability criterion
	if ρ_ex == Inf
		circles = gershgorin!(linearize(y -> i.f_ex(i.t, y), i.y), (i.w_1, i.w_2, i.w_3))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end

	i.dt_[] = safety * (3/sqrt2) / ρ_ex
	i.dt[]  = min(i.dt[], i.dt_[])
	
	assign!(i.w_1, i.y)

	@algo_step begin
		@verbo_println("first trapezoidal stage:")
		assign!(i.w_2, i.tr1_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_2)
			sc = SchurComplement(w -> i.rhs - i.tr1_lhs(w), i.w_2[2])
			newtonit!(
				sc, i.w_2[1], i.h[1], i.h[2:end];
				maxit    = newton_maxit,
				rtol     = newton_rtol,
				cg_maxit = cg_maxit,
				cg_rtol  = cg_rtol
			)
		end
	end

    @algo_step begin
		@verbo_println("bdf2 stage:")
		assign!(i.w_3, i.bdf2_rhs)
    	if i.imex.im
			assign!(i.rhs, i.w_3)
			sc = SchurComplement(w -> i.rhs - i.bdf2_lhs(w), i.w_3[2])
			newtonit!(
				sc, i.w_3[1], i.h[1], i.h[2:end];
				maxit    = newton_maxit,
				rtol     = newton_rtol,
				cg_maxit = cg_maxit,
				cg_rtol  = cg_rtol
			)
		end
	end

	@algo_step begin
		@verbo_println("second trapezoidal stage:")
		assign!(i.w_4, i.tr2_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_4)
			sc = SchurComplement(w -> i.rhs - i.tr2_lhs(w), i.w_4[2])
			newtonit!(
				sc, i.w_4[1], i.h[1], i.h[2:end];
				maxit    = newton_maxit,
				rtol     = newton_rtol,
				cg_maxit = cg_maxit,
				cg_rtol  = cg_rtol
			)
		end
	end

	assign!(i.e, i.err_ex + i.err_im)
	
	η = l2.(i.e) / (rtol * l2.(i.y) + atol)

	@verbo_println("TR-BDF2-SCR $(i.k[]), t = $(@sprintf "%.3e" i.t[]), dt = $(@sprintf "%.3e" i.dt[]), η = $(@sprintf "(%.3e %.3e)" η...)")
	
	if i.k[] == 0 && i.dt[] < i.dt_[] && abs(max(η...) - 1) > .1
		i.dt[] *= 0.1 + 0.9/sqrt(max(η...))
		return step!(i;
			rtol         = rtol,
			atol         = atol,
			newton_maxit = newton_maxit,
			newton_rtol  = newton_rtol,
			cg_maxit     = cg_maxit,
			cg_rtol      = cg_rtol,
			growth       = growth,
			safety       = safety,
			ρ_ex         = ρ_ex
		)
	end

	assign!(i.y, i.update_y)
	
	i.t[] += i.dt[]
	i.k[] += 1

	i.dt[] = min(1/sqrt(max(η...)), growth) * i.dt[]

	return return (; zip((:v, :e), η)...)
end

function visualize_op(op, vecs)
	assign!(vecs[2], diag(vecs[1], op(vecs[1]), (-1,)))
	assign!(vecs[3], diag(vecs[1], op(vecs[1]), ( 0,)))
	assign!(vecs[4], diag(vecs[1], op(vecs[1]), (+1,)))
	assign!(vecs[2], vecs[2]/vecs[3])
	assign!(vecs[4], vecs[4]/vecs[3])
	Meta.@show vecs[2].y.data[2:end]
	Meta.@show vecs[4].y.data[1:end-1]
	display <| plot([vecs[2].y.data[2:end] vecs[4].y.data[1:end-1]]; label=["-1" "+1"])
end

# Can also be used for breakage only
function step!(i::tr_bdf2{Intg_DR}; rtol, atol = 0, newton_maxit, newton_rtol, cg_maxit, cg_rtol, growth, safety, ρ_ex = Inf)
	
	if ρ_ex == 0 && i.dt[] == Inf # fully implicit
		i.dt[] = 1 # some value, to be calibrated
	end

	if ρ_ex == Inf
		circles = gershgorin!(linearize(((e, α),) -> i.f_ex(i.t[], (i.y[1], i.y[2], e, α),)[3:4], i.y[3:4]), (i.w_1[3:4], i.w_2[3:4], i.w_3[3:4]))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end
	
	# TODO: improve stability criterion
	i.dt_[] = safety * (3/sqrt2)  / ρ_ex
	i.dt[]  = min(i.dt[], i.dt_[])
	
	assign!(i.w_1, i.y)

	@algo_step begin
		@verbo_println("first trapezoidal stage:")
		assign!(i.w_2, i.tr1_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_2)
			
			@algo_step begin
				@verbo_println("elasticity:")
				elastic = SchurComplement(((v, e),) -> i.rhs[2:3] - i.tr1_lhs((i.y[1], v, e, i.y[4]),)[2:3], i.w_2[3])
				
				# visualize_op(linearize(elastic, i.w_2[1]), i.h.v)
				# error("asdf")
				
				newtonit!(
					elastic, i.w_2[2], i.h.v[1], i.h.v[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
			
			@algo_step begin
				@verbo_println("damage or breakage:")
				damage  = α -> i.rhs[4] - i.tr1_lhs((i.y[1], i.y[2], i.y[3], α),)[4]
				newtonit!(
					damage, i.w_2[4], i.h.α[1], i.h.α[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
		end
	end

	i.t[] += (2 - sqrt2) * i.dt[]

	@algo_step begin
		@verbo_println("bdf2 stage:")
    	assign!(i.w_3, i.bdf2_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_3)
			
			@algo_step begin
				@verbo_println("elasticity:")
				elastic = SchurComplement(((v, e),) -> i.rhs[2:3] - i.bdf2_lhs((i.y[1], v, e, i.y[4]),)[2:3], i.w_3[3])
				newtonit!(
					elastic, i.w_3[2], i.h.v[1], i.h.v[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
			
			@algo_step begin
				@verbo_println("damage or breakage:")
				damage  = α -> i.rhs[4] - i.bdf2_lhs((i.y[1], i.y[2], i.y[3], α),)[4]
				newtonit!(
					damage,  i.w_3[4], i.h.α[1], i.h.α[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
		end
	end

	i.t[] -= (2 - sqrt2) * i.dt[]
	
	@algo_step begin
		@verbo_println("second trapezoidal stage:")
		assign!(i.w_4, i.tr2_rhs)
		if i.imex.im
			assign!(i.rhs,  i.w_4)
			
			@algo_step begin
				@verbo_println("elasticity:")
				elastic = SchurComplement(((v, e),) -> i.rhs[2:3] - i.tr2_lhs((i.y[1], v, e, i.y[4]),)[2:3], i.w_4[3])
				newtonit!(
					elastic, i.w_4[2], i.h.v[1], i.h.v[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
			
			@algo_step begin
				@verbo_println("damage or breakage:")
				damage  = α -> i.rhs[4] - i.tr2_lhs((i.y[1], i.y[2], i.y[3], α),)[4]
				newtonit!(
					damage,  i.w_4[4], i.h.α[1], i.h.α[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
		end
	end

	assign!(i.e[1], i.err_im[1] + i.err_ex[1])
	assign!(i.e[2], i.err_im[2] + i.err_ex[2])
	assign!(i.e[3], i.err_im[3] + i.err_ex[3])
	assign!(i.e[4], i.err_im[4] + i.err_ex[4])
	
	# note: relative to w_3 instead of new y
	η = l2.(i.e) ./ (values(rtol) .* l2.(i.w_3) .+ values(atol))

	@verbo_println("TR-BDF2-DR $(i.k[]), t = $(@sprintf "%.3e" i.t[]), dt = $(@sprintf "%.3e" i.dt[]), η = $(@sprintf "(u = %.3e, v = %.3e, e = %.3e, α = %.3e)" η[1] η[2] η[3:4]...)")
	
	if i.k[] == 0 && i.dt[] < i.dt_[] && abs(max(values(η)...) - 1) > .1
		i.dt[] *= 0.1 + 0.9/sqrt(max(values(η)...))
		return step!(i;
			rtol         = rtol,
			atol         = atol,
			newton_maxit = newton_maxit,
			newton_rtol  = newton_rtol,
			cg_maxit     = cg_maxit,
			cg_rtol      = cg_rtol,
			growth       = growth,
			safety       = safety,
			ρ_ex         = ρ_ex
		)
	end

	assign!(i.y, i.update_y)
	
	i.t[] += i.dt[]
	i.k[] += 1

	i.dt[] = min(1/sqrt(max(η...)), growth) * i.dt[]

	return (; zip((:u, :v, :e, :α), η)...)
end


function step!(i::tr_bdf2{Intg_DBR}; rtol, atol = 0, newton_maxit, newton_rtol, cg_maxit, cg_rtol, growth, safety, ρ_ex = Inf)
	
	if ρ_ex == 0 && i.dt[] == Inf # fully implicit
		i.dt[] = 1 # some value, to be calibrated
	end

	if ρ_ex == Inf
		circles = gershgorin!(linearize(((e, α, β),) -> i.f_ex(i.t[], (i.y[1], e, α, β),)[2:4], i.y[2:4]), (i.w_1[2:4], i.w_2[2:4], i.w_3[2:4]))
		ρ_ex = max(map(c -> abs(c.o .- c.r), circles)...)
	end

	# TODO: improve stability criterion
	i.dt_[] = safety * (3/sqrt2) / ρ_ex
	i.dt[]  = min(i.dt[], i.dt_[])
	
	assign!(i.w_1, i.y)

	@algo_step begin
		@verbo_println("first trapezoidal stage:")
		assign!(i.w_2, i.tr1_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_2)
			
			@algo_step begin
				@verbo_println("elasticity:")
				elastic  = SchurComplement(((v, e),) -> i.rhs[1:2] - i.tr1_lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_2.e)
				newtonit!(
					elastic, i.w_2.v, i.h.v[1], i.h.v[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
			
			@algo_step begin
				@verbo_println("damage:")
				damage   = α -> i.rhs[3] - i.tr1_lhs((i.y.v, i.y.e, α, i.y.β),)[3]
				newtonit!(
					damage,   i.w_2.α, i.h.α[1], i.h.α[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end

			@algo_step begin
				@verbo_println("breakage:")
				breakage = β -> i.rhs[4] - i.tr1_lhs((i.y.v, i.y.e, i.y.α, β),)[4]
				newtonit!(
					breakage, i.w_2.β, i.h.β[1], i.h.β[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
		end
	end

	@algo_step begin
		@verbo_println("bdf2 stage:")
    	assign!(i.w_3, i.bdf2_rhs)
		if i.imex.im
			assign!(i.rhs, i.w_3)
			
			@algo_step begin
				@verbo_println("elasticity:")
				elastic = SchurComplement(((v, e),) -> i.rhs[1:2] - i.bdf2_lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_3.e)
				newtonit!(
					elastic,  i.w_3.v, i.h.v[1], i.h.v[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
			
			@algo_step begin
				@verbo_println("damage:")
				damage   = α -> i.rhs[3] - i.bdf2_lhs((i.y.v, i.y.e, α, i.y.β),)[3]
				newtonit!(
					damage,   i.w_3.α, i.h.α[1], i.h.α[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
			
			@algo_step begin
				@verbo_println("breakage:")
				breakage = β -> i.rhs[4] - i.bdf2_lhs((i.y.v, i.y.e, i.y.α, β),)[4]
				newtonit!(
					breakage, i.w_3.β, i.h.β[1], i.h.β[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
		end
	end
	
	@algo_step begin
		@verbo_println("second trapezoidal stage:")
		assign!(i.w_4, i.tr2_rhs)
		if i.imex.im
			assign!(i.rhs,  i.w_4)
			
			@algo_step begin
				@verbo_println("elasticity:")
				elastic  = SchurComplement(((v, e),) -> i.rhs[1:2] - i.tr2_lhs((v, e, i.y.α, i.y.β),)[1:2], i.w_4.e)
				newtonit!(
					elastic,  i.w_4.v, i.h.v[1], i.h.v[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end

			@algo_step begin
				@verbo_println("damage:")
				damage   = α -> i.rhs[3] - i.tr2_lhs((i.y.v, i.y.e, α, i.y.β),)[3]
				newtonit!(
					damage,   i.w_4.α, i.h.α[1], i.h.α[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end

			@algo_step begin
				@verbo_println("breakage:")
				breakage = β -> i.rhs[4] - i.tr2_lhs((i.y.v, i.y.e, i.y.α, β),)[4]
				newtonit!(
					breakage, i.w_4.β, i.h.β[1], i.h.β[2:end];
					maxit    = newton_maxit,
					rtol     = newton_rtol,
					cg_maxit = cg_maxit,
					cg_rtol  = cg_rtol
				)
			end
		end
	end

	assign!(i.e, i.update_e)
	
	# note: relative to w_3 instead of new y
	η = l2.(i.e) ./ (rtol .* l2.(i.w_3) .+ atol)

	@verbo_println("TR-BDF2-DBR $(i.k[]), t = $(i.t[]), dt = $(i.dt[]), η = $η")
	
	if i.k[] == 0 && i.dt[] < i.dt_[] && abs(max(values(η)[2:4]...) - 1) > .1
		i.dt[] *= 0.1 + 0.9/sqrt(max(values(η)[2:4]...))
		return step!(i;
			rtol         = rtol,
			atol         = atol,
			newton_maxit = newton_maxit,
			newton_rtol  = newton_rtol,
			cg_maxit     = cg_maxit,
			cg_rtol      = cg_rtol,
			growth       = growth,
			safety       = safety,
			ρ_ex         = ρ_ex
		)
	end

	assign!(i.y, i.update_y)
	
	i.t[] += i.dt[]
	i.k[] += 1

	i.dt[] = min(1/sqrt(max(η...)), growth) * i.dt[]

	return (; zip((:v, :e, :α, :β), η)...)
end

