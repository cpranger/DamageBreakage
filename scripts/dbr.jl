#!/usr/bin/env -S julia --project

name = "dbr-plastic"

const Float = Float32

const JULIA_CUDA_USE_BINARYBUILDER = false
const USE_GPU = false
const GPU_ID  = 1

# using ParallelStencil
# # using ParallelStencil.FiniteDifferences2D
# @static if USE_GPU
# 	@init_parallel_stencil(CUDA, Float, 2)
# 	CUDA.device!(GPU_ID)
# else
# 	@init_parallel_stencil(Threads, Float, 2)
# end
using JLD
using Printf
# using Adapt

using StaggeredKernels
using StaggeredKernels.Plane

include("../algorithm/algorithm.jl")

using .Algo

Strain(dims)   = StaggeredKernels.Tensor(dims, StaggeredKernels.Symmetric,      strain_stags...)
Velocity(dims) = StaggeredKernels.Tensor(dims, StaggeredKernels.Unsymmetric{1}, motion_stags...)
Damage(dims)   = StaggeredKernels.Field(dims, state_stags...)

# aliases
const Breakage = Damage
const Norm     = Damage
const Coord    = Velocity

State(n) = (v = Velocity(n), e = Strain(n), α = Damage(n), β = Breakage(n))

# Utils

softramp(ρ, x) = If(real(x) >= ρ, x, ρ*exp(x/ρ-1))

grad_sq(a) = Tensor((
	xx = A(D(a, :x), :y) * A(D(a, :x), :y),
	xy = A(D(a, :x), :y) * A(D(a, :y), :x),
	yy = A(D(a, :y), :x) * A(D(a, :y), :x),
), Symmetric)

diag_lap(f) = 2*(4*A(A(f, :y), :x) - 4*f)


# Physics

ξ(e)    = I1(e) / sqrt(I2(e))
f(e, p) = sqrt(J2(e) / softramp(p.ε_0^2, max(0, real(p.e_0 - I1(e)))^2 - J2(e))) / p.f_0

λ(e, α, p) = p.λ_0 - α * p.γ_r / ξ(e)

# note: factor 2 added in last term w.r.t. Lyakhovksy and Ben-Zion, 2014b.
μ(e, α, p) = p.μ_0 - α * p.γ_r * (ξ(e) - p.ξ_0)

C(e, α, p) = λ(e, α, p) * MajorIdentity() + μ(e, α, p) * MinorIdentity()

s(e, α, p) = C(e, α, p) * e - p.θ_0 * grad_sq(α) / p.h^2

γ(e, β, p) = f(e, p)^p.f_a * β^p.f_b

# only _elastic_ density changes
r(e, p) = p.r_0 * exp(-I1(e))

# Schmidt tensor
S(s) = dev(s) / sqrt(J2(s))
S(e, α, p) = S(s(e, α, p))

F_im(v, e, α, β, p) = (
    v =   divergence(s(e, α, p)) / p.h / interpolate(r(e, p)),
    e =   symgrad(v) / p.h,
    α =   0,  # TODO
    β = - #=scale *=# U.β * p.l_0^2 * diag_lap(γ(e, β, p)) / p.h^2,
)

F_ex(v, e, α, β, p) = (
    v =   (x = 0, y = 0),
    e = - p.γ_0 * γ(e, β, p) * S(e, α, p),
    α =   0,  # TODO
    β = - #=scale *=# β * γ(e, β, p)
)

F(U, p) = F_ex(U, p) + F_im(U, p)


function random_init()
	gen = StaggeredKernels.FieldGen((i, j) -> rand())
	return (
	    v = (x = gen, y = gen),
	    e = (
	        xx = gen,
	        xy = gen,
	        yy = gen
	    ),
	    α = gen,
	    β = gen
	)
end


compute_params(p) = let λ_0 = p.μ_0   # assuming Poisson's ratio ν = 1/4
	merge(p, (
	l_x   = (p.n_x-2)*p.h,            # domain size (1)
	l_y   = (p.n_y-2)*p.h,            # domain size (2)
	γ_0   = p.v_0/p.l_0,                # reference strain rate (/s)
	λ_0   = λ_0,                        # Lamé's first parameter 
	k_0   = λ_0 + (2/3)*p.μ_0,          # Bulk modulus
	c_s   = sqrt(       p.μ_0 /p.r_0),  # solenoidal  wave speed (m/s)
	c_p   = sqrt((λ_0+2*p.μ_0)/p.r_0),  # compressive wave speed (m/s)
	range = (1:p.n_x, 1:p.n_y)
)) end



λ(i::Int,    n, b, s, h) = 2/h * sin(h/2 * ν(i, n, b, s, h))
ν(i::Int,    n, b, s, h) = π/h * (1 - 1/2*(b['-'] - b['+'])^2) /
                                 (n - 1/2 - (1 - b['-'] - b['+'])*(1/2 - s))
φ(i::Int, j, n, b, s, h) = sin(ν(i, n, b, s, h) * h * (j + (1/2 - s)*b['-']) + π/2 * b['-'])

κ(k, i::Vector, n, b_l, b_u, s, h) = (-1)^b[k][k]['-'] * λ(i[k], n[k], b[k][k], s[k][k], h)

φ(i::Vector, j, n, b, s, h) = *(
	φ(i[1], j[1], n[1], b[1], s[1], h),
	φ(i[2], j[2], n[2], b[2], s[2], h),
	φ(i[3], j[3], n[3], b[3], s[3], h)
)

φ_p(i::Vector, j, n, b, s, h) = [
	κ(1, i, n, b, s, h) * φ(i, n, b[1], s[1], h),
	κ(2, i, n, b, s, h) * φ(i, n, b[2], s[2], h),
	κ(3, i, n, b, s, h) * φ(i, n, b[3], s[3], h)
]

φ_s(i::Vector, j, n, b, s, h) = [
	(- κ(3, i, n, b, s, h) + κ(2, i, n, b, s)) * φ(i, n, b[1], s[1], h),
	(+ κ(3, i, n, b, s, h) - κ(1, i, n, b, s)) * φ(i, n, b[2], s[2], h),
	(- κ(2, i, n, b, s, h) + κ(1, i, n, b, s)) * φ(i, n, b[3], s[3], h)
]

join_b(b_l, b_u) = map(
	(b_l, b_u) ->  map(
		(b_l, b_u) ->  Dict(['-' => b_l, '+' => b_u]),
		 b_l, b_u
	 ),
	 b_l, b_u
)

function eigenmode_v(p)
	i = [2, 3, 0]
	# rows: dimensions
	# cols: components
	s = [
		[0, 1, 1],  # dimension 1
		[1, 0, 1],  # dimension 2
		[1, 1, 0],  # dimension 3
	]
	b_l = [
		[0, 1, 0],  # dimension 1
		[1, 0, 0],  # dimension 2
		[1, 1, 1],  # dimension 3
	]
	b_u = [
		[0, 1, 0],  # dimension 1
		[1, 0, 0],  # dimension 2
		[1, 1, 1],  # dimension 3
	]
	n = [
		p.n_x + 1 .- s[1],
		p.n_y + 1 .- s[2],
		p.n_z + 1 .- s[3],
	]
	h = p.h
	
	return (
		fieldgen((j1, j2) -> φ_p(i, [j1-1, j2-1, 0], n, b, h)),
		fieldgen((j1, j2) -> φ_s(i, [j1-1, j2-1, 0], n, b, h))
	)
end

function dbr()
	# TODO: get these from options!
	blocksize = 10 # 64    # size of an array block
	nblocks = (1, 1)       # number of blocks in each dimension
	dims  = (n_x, n_y) = nblocks .* blocksize  # resolution
	
	yr = Float(π * 1.e7)    # year (s)
	
	l_x        =   10.e3
	h          =   l_x/(n_x-2)
	
	p = compute_params((
		n_x    =   n_x   ,               # resolution (1)
		n_y    =   n_y   ,               # resolution (2)
		h      =   h     ,               # cell size (isometric)
		f_0    =   0.6   ,               # ref. friction coeff. (-)
		f_a    =   0.001 ,               # a strain rate exponent (-)
		f_b    =   0.001 ,               # b grain  size exponent (-)
		l_0    =   2*h ,                 # diffusion length scale (m)
		v_0    =   1.e-9 ,               # reference velocity (m/s)
		d_0    =   1.e-3 ,               # reference slip (m)
		μ_0    =   5e9   ,               # Lamé's second parameter | shear modulus (Pa)
		r_0    =   2500. ,               # density (kg/m^3)
		e_0    =   5.e-4 ,               # cohesive volumetric strain (-)
		ε_0    =   1.e-4 ,               # regularization strain
		γ_r    =   0.    ,               # elastic damage coefficient 1
		ξ_0    =   0.    ,               # elastic damage coefficient 2
		θ_0    =   0.    ,               # elastic damage coefficient 3
		h_t    =   0.                    # time step
	))
	
	# Fields
	X         =  Coord(dims)  # coordinate
	U         =  State(dims)  # state vector
	Λ         =  State(dims)  # eigenvector
	H         = [State(dims) for _ in 1:6] # Helper vectors
	H2        =   Norm(dims)               # Helper vector
	
	assign!(
	    X, (
		     x = fieldgen((i,j) -> i * p.h - p.l_x/2),
		     y = fieldgen((i,j) -> j * p.h - p.l_y/2)
		), p.range
	)
	
	# initial conditions
	assign!(U, (
	    v = eigenmode_v(),
	    e = eigenmode_e(),
	    α = 0,
	    β = 0, # If(L2(X) < 2*π*p.l_0, #=(scale TODO) * =# cos(0.25/p.l_0*L2(X))^2)
	), p.range)
	
	
	# randomly initialize Λ
	assign!(Λ, random_init(), p.range)
	
	assign!(Λ, F_ex(U, p), p.range)
	
	assign!(Λ, F_im(U, p), p.range)
	
	# step = 0; nstep = 1# ; fr = 0
	# while step < nstep
	#
	# 	# TODO: select explicit parts of U
	# 	# 1 compute eigenvalue of explicit part of F
	# 	eigen!(Λ, U, (H2, H[1], H[2]), p, 1e-6 #=rtol=#)
	#
	# 	# 2 compute time step from that
	# 	p.h_t = sqrt(3) / max(H2) # TODO: Check!
	#
	# 	# 3 make time step
	# 	trbdf2_step!(U, (H[1], H[2], H[3]), H[4], H[5], H[6], p)
	#
	# 	# 4 increment counter
	# 	step += 1
	# end
	
	"finished!"
end

dbr()