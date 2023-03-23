# module Algo

cstep_sow( a, b, p) = a + p.ħ * im * b
cstep_reap(a, p) = Im(a) / p.ħ

# TR-BDF2 SCHEME

# LHS
trbdf2_lhs(W_n, p) = W_n - p.h_t * ((1-sqrt(2)/2) * F_im(W_n, p))

# TR-BDF2 stage 1,2,E RHS
trbdf2_rhs(Y, W_0, W_1, W_2, p) = (
    Y + p.h_t * (
        (2/1-sqrt(2)/1) * F_ex(W_0, p)
      + (1/1-sqrt(2)/2) * F_im(W_0, p)
	),
    Y + p.h_t * (
        (1/2-sqrt(2)/3) * F_ex(W_0, p) + (1/2+sqrt(2)/3) * F_ex(W_1, p)
      + (    sqrt(2)/4) * F_im(W_0, p) + (    sqrt(2)/4) * F_im(W_1, p)
    ),
	    p.h_t * (
	    (1/3-sqrt(2)/3) * F_im(W_0, p) + (1/3) * F_im(W_1, p) + (-2/3 + sqrt(2)/3)*F_im(W_2, p)
	)
)

# TR-BDF2 update
trbdf2_upd(Y, W_0, W_1, W_2, p) = Y + p.h_t * (
    (sqrt(2)/4)*F(W_0, p) + (sqrt(2)/4)*F(W_1, p) + (1/1-sqrt(2)/2)*F(W_2, p)
)

trbdf2_implicit_JV(U, V, p) = (
    v = (
        x = cstep_reap(trbdf2_lhs((v = U.v, e = cstep_sow(U.e, V.e), α = U.α, β = U.β), p).v.x),
        y = cstep_reap(trbdf2_lhs((v = U.v, e = cstep_sow(U.e, V.e), α = U.α, β = U.β), p).v.y)
    ),
    e = (
        xx = cstep_reap(trbdf2_lhs((v = cstep_sow(U.v, V.v), e = U.e, α = U.α, β = U.β), p).e.xx),
        xy = cstep_reap(trbdf2_lhs((v = cstep_sow(U.v, V.v), e = U.e, α = U.α, β = U.β), p).e.xy),
        yy = cstep_reap(trbdf2_lhs((v = cstep_sow(U.v, V.v), e = U.e, α = U.α, β = U.β), p).e.yy)
    ),
    α = cstep_reap(trbdf2_lhs((v = U.v, e = U.e, α = cstep_sow(U.α, V.α), β = U.β), p).α),
    β = cstep_reap(trbdf2_lhs((v = U.v, e = U.e, α = U.α, β = cstep_sow(U.β, v.β)), p).β)
)

function trbdf2_residual_v!(R, U, V, B, p)
	# schur decomposition of the (v, e)-subproblem
	assign!(V.e,       trbdf2_implicit_JV(U, V, p).e,       p.ranges)
	assign!(R.v, V.v - trbdf2_implicit_JV(U, V, p).v - B.v, p.ranges)
end

trbdf2_residual_α!(R, U, V, B, p) = assign!(R.α, trbdf2_implicit_JV(U, V, p).α - B.α, p.ranges)
trbdf2_residual_β!(R, U, V, B, p) = assign!(R.β, trbdf2_implicit_JV(U, V, p).β - B.β, p.ranges)

function trbdf2_stage!(R, Y, U, V, B, RHS, p)
	assign!(B, RHS, p.ranges)
	
	# schur decomposition: b_v = b_v - B b_e
	assign!(B.v, B.v - trbdf2_implicit_JV(U, B, p).v, p.ranges)
	
	Chebyshev!((lhs, rhs) -> assign!(lhs.v, rhs.v), trbdf2_residual_v!, R, U, V, B, #=...,=# p, p.k_1_v, p.k_n_v)
	Chebyshev!((lhs, rhs) -> assign!(lhs.α, rhs.α), trbdf2_residual_α!, R, U, V, B, #=...,=# p, p.k_1_α, p.k_n_α)
	Chebyshev!((lhs, rhs) -> assign!(lhs.β, rhs.β), trbdf2_residual_β!, R, U, V, B, #=...,=# p, p.k_1_β, p.k_n_β)
	
	# schur decomposition: v_e = b_e - C v_v
	assign!(V.e, B.e - trbdf2_implicit_JV(U, V, p).e, p.ranges)
end

trbdf2_finish!(Y, W_0, W_1, W_2, p) = assign!(Y, trbdf2_upd(W_0, W_1, W_2, p), p.ranges)

function trbdf2_step!(Y, E, W_0, W_1, W_2, B, R, p)
	RHS = trbdf2_rhs(Y, W_0, W_1, W_2, p)
	
	# stage 0
	assign!(W_0, Y, p)
	
	# stage 1
	assign!(W_1, W_0, p)  # initial guess
	trbdf2_stage!(R, Y, W_0, W_1, B, RHS[1], p)
	
	# stage 2
	assign!(W_2, W_1, p)  # initial guess
	trbdf2_stage!(R, Y, W_1, W_2, B, RHS[2], p)
	
	# error stage
	trbdf2_stage!(R, Y, W_2, E,   B, RHS[3], p)
	
	# finish
	trbdf2_finish!(Y, W_0, W_1, W_2, p)
end


# function Chebyshev!(assign, residual!, R, U, X, B, V, p, k_1, k_n)
	
# 	# based on Gutknecht & Röllin (2002; Parallel Computing), algorithm 5.
# 	Err  = []
	
# 	# r = [J]_u * x - b
# 	residual(R, U, X, B, p)
	
# 	push!(Err, norm(r.c) / k_1)
	
# 	if Err[1] < atol return Err end
	
# 	α = -(k_n + k_1) / 2
# 	c = -(k_n - k_1) / 2
	
# 	cheby = (i,x::Complex) -> real(0.5*((x + sqrt(x^2 - 1))^i + (x - sqrt(x^2 - 1))^i))
	
# 	if abs(c/α) < atol/Err[1]
# 		ρ = i -> abs(c/α)
# 	else
# 		ρ = i -> 1/cheby(i,α/c + 0im)
# 	end
	
# 	# massive overstimate, as if simple Jacobi iteration
# 	nit_im = max(2,ceil(UInt64, log(atol/Err[1])/log(ρ(1))))
	
# 	ψ = 0
# 	ω = 1/α
	
# 	i = 0; while i < nit_im
# 		if i == 1
# 			ψ =      + (1/2)*(c/α)^2  # - -> +
# 			ω = 1/(α - c^2/(2*α))
# 		elseif i > 1
# 			ψ =      + (c/2)^2 * ω^2  # - -> +
# 			ω = 1/(α - (c/2)^2 * ω)
# 		end
		
# 		# TODO: allow arithmetics on objects of type Unknown
# 		assign!(V, R + ψ*V, p.range)  # - -> +
# 		assign!(X, X + ω*V, p.range)
		
# 		residual!(R, U, X, B, p)
		
# 		i += 1
		
# 		# TODO perf: do only sparsly
# 		# push!(Err, sqrt(norm(r.x)^2 + norm(r.y)^2) / k_1)
# 		push!(Err, Err[1]*ρ(i))
		
# 		# if Err[1]*ρ(i) < atol break end
# 		if Err[end] < atol break end
# 	end
	
# 	return Err
# end

# rayleigh_quotient!(Λ, JB, B, p) = assign!(Λ, dot(JB, B) / dot(B, B), p.ranges)

# function eigen!(B, U, Λ, E, p)
# 	for _ in 1:100
# 		JB = cstep_reap(F_ex(cstep_sow(U, B, p), p), p)
		
# 		# λ = Jb.b / b.b
# 		local_rayleigh_quotient!(Λ, JB, B, p)
		
# 		# err = || Jb / b - λ ||
# 		ErrExpr = (
# 		    e = (
# 		        xx = JB.e.xx / B.e.xx - Λ,
# 		        xy = JB.e.xy / B.e.xy - Λ,
# 		        yy = JB.e.yy / B.e.yy - Λ
# 		    ),
# 			α = JB.α / B.α - Λ,
# 		    β = JB.β / B.β - Λ
# 		)
# 		assign_expl!(E, sqrt(local_dot_expl(ErrExpr, ErrExpr)))
# 		if norm(E)/norm(Λ) < p.rtol break end
		
# 		# renormalize b
# 		assign_expl!(B, B/sqrt(local_dot_expl(B, B)))
# 	end
# end

# end # module Algo