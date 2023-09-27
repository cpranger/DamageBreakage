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
	assign!(V.e,       trbdf2_implicit_JV(U, V, p).e)
	assign!(R.v, V.v - trbdf2_implicit_JV(U, V, p).v - B.v)
end

trbdf2_residual_α!(R, U, V, B, p) = assign!(R.α, trbdf2_implicit_JV(U, V, p).α - B.α)
trbdf2_residual_β!(R, U, V, B, p) = assign!(R.β, trbdf2_implicit_JV(U, V, p).β - B.β)

function trbdf2_stage!(R, Y, U, V, B, RHS, p)
	assign!(B, RHS)
	
	# schur decomposition: b_v = b_v - B b_e
	assign!(B.v, B.v - trbdf2_implicit_JV(U, B, p).v)
	
	Chebyshev!((lhs, rhs) -> assign!(lhs.v, rhs.v), trbdf2_residual_v!, R, U, V, B, #=...,=# p, p.k_1_v, p.k_n_v)
	Chebyshev!((lhs, rhs) -> assign!(lhs.α, rhs.α), trbdf2_residual_α!, R, U, V, B, #=...,=# p, p.k_1_α, p.k_n_α)
	Chebyshev!((lhs, rhs) -> assign!(lhs.β, rhs.β), trbdf2_residual_β!, R, U, V, B, #=...,=# p, p.k_1_β, p.k_n_β)
	
	# schur decomposition: v_e = b_e - C v_v
	assign!(V.e, B.e - trbdf2_implicit_JV(U, V, p).e)
end

trbdf2_finish!(Y, W_0, W_1, W_2, p) = assign!(Y, trbdf2_upd(W_0, W_1, W_2, p))

function trbdf2_step!(Y, E, W_0, W_1, W_2, B, R, p)
	RHS = trbdf2_rhs(Y, W_0, W_1, W_2, p)
	
	# stage 0
	assign!(W_0, Y)
	
	# stage 1
	assign!(W_1, W_0)  # initial guess
	trbdf2_stage!(R, Y, W_0, W_1, B, RHS[1], p)
	
	# stage 2
	assign!(W_2, W_1)  # initial guess
	trbdf2_stage!(R, Y, W_1, W_2, B, RHS[2], p)
	
	# error stage
	trbdf2_stage!(R, Y, W_2, E,   B, RHS[3], p)
	
	# finish
	trbdf2_finish!(Y, W_0, W_1, W_2, p)
end
