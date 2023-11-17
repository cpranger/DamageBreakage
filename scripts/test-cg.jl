#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane
# using LinearAlgebra

function test_poisson(p, ax)
	u = Field(p.n, div_stags)
	h = Tuple([deepcopy(u) for _ in 1:5])
	b = deepcopy(u)
	
	f  = u -> divergence(grad(u))# - 3*p.h^2
	bc(u) = (
		"-y" =>  FD(u, :y),
		"+y" => 1 - u,
		"-x" =>    -u,
		"+x" =>    -u,
	)

	A = linearize(u -> (f(u), bc(u)), 0)
	assign!(b, (f(0), bc(0)))
	assign!(b, -b)
	
	r = deepcopy(u)
	
	# Meta.@show n = p.n[1]

	# B = zeros(n, n, n, n)

	# assign!(u, 0)
	# for i in 1:n
	# 	for j in 1:n
	# 		for k in 1:n
	# 			for l in 1:n
	# 				u.data[1, k, l] = 1;
	# 				assign!(r, A(u))
	# 				B[i, j, k, l] = r.data[1, i, j]
	# 				u.data[1, k, l] = 0;
	# 			end
	# 		end
	# 	end
	# end

	# C = reshape(B, (n^2, n^2))
	# display(C)

	# d = reshape(b.data, n^2)
	# display(d)

	# display <| Plots.heatmap(1:n^2, 1:n^2, C, aspectratio = :equal, framestyle  = :box)

	# println("Eigenvalues")
	# display <| LinearAlgebra.eigvals(C)
	# println("Eigenvectors")
	# display <| LinearAlgebra.eigvecs(C)

	# Meta.@show C\d

	# Meta.@show LinearAlgebra.isposdef(.-C)

	# # return
	# readline()

	assign!(u, (0, bc(0)))
	# display <| reshape(u.data, n^2)

	# r = deepcopy(u)
	# assign!(r, b - A(u))
	# display <| r.data
	
	# (λ, Λ, ε) = cg!(A, u, b; h = h[1:3], rtol = 1e-6, λtol = 1e-2, minit = 100, maxit = 1000)
	(λ, Λ, ε) = cg_pc_jacobi!(A, u, b; h = h, rtol = 1e-6, λtol = 1e-2, minit = 12, maxit = 100)
	# println  <| cg_petsc!(A, u, b, h, ax; normtype = :norm_preconditioned, hermitian = false, rtol = 1e-6, atol = 1e-6, dtol = 1e1, maxit = 100)
	
	r = deepcopy(u)
	assign!(r, b - A(u))

	plt1 = heatmap(ax..., u, "u", c = :davos)
	plt2 = heatmap(ax..., r, "r", c = :davos)
	display <| plot(plt1, plt2; layout = (1, 2))
	
	return
end

function test_elastic(p, ax_x, ax_y)
	u   = Vector(p.n, motion_stags)
	b   = deepcopy(u)
	h   = Tuple([deepcopy(u) for _ in 1:5])
	
	s(e) = (p.λ_0 * Tensor(MajorIdentity) + p.μ_0 * Tensor(MinorIdentity)) * e
	
	f  = u -> divergence(s(symgrad(u)))
	bc = u -> (;
		x = (
			"-y" => FD(  u.x, :y),
			"-x" =>   ( -u.x    ),
			"+x" =>   ( -u.x    ),
			"+y" =>   (fieldgen((i, j) -> sin(pi*i/(p.n[1]-2))) - u.x)
		),
		y = (
			"-y" =>   ( -u.y    ),
			"-x" =>   ( -u.y    ),
			"+x" =>   ( -u.y    ),
			"+y" =>   ( -u.y    )
		)
	)
	
	A = linearize(u -> (f(u), bc(u)), 0*u)
	assign!(b, (f(u), bc(u)))
	assign!(b, -b)
	
	assign!(u, (0, bc(0*u)))
	(λ, Λ, ε) = cg!(A, u, b; h = h[1:3], rtol = 1e-8, λtol = 1e-2, minit = 100, maxit = 1000)
	# (λ, Λ, ε) = cg_pc_jacobi!(A, u, b; h = h, rtol = 1e-8, λtol = 1e-2, minit = 100, maxit = 1000)

	r = h[2]; assign!(r, b - A(u))

	plt1 = heatmap(ax_x, ax_y, u, "u", c = :davos)
	plt2 = heatmap(ax_x, ax_y, r, "r", c = :davos)
	display <| plot(plt1, plt2; layout = (1, 2))
	
	readline()
	return
end

function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	h    =  1 / (n[1] - 2)
	
	μ_0  =  1
	λ_0  =  μ_0                      # assuming Poisson's ratio ν = 1/4
	k_0  =  λ_0 + (2/3)*μ_0          # bulk modulus

	# collect all variables local to this function:
	vars = Base.@locals
	
	# and return as a named tuple
	return NamedTuple{Tuple(keys(vars))}(values(vars))
end

function main()
	s = ArgParseSettings()
	
	@add_arg_table s begin
    	"--nb"
			help = "comma-separated list of the number of blocks in each dimension"
			required = true
			arg_type = (NTuple{N, Int} where N)
    end
	
	p = parameters(; parse_args(s)...)

	# axes
	ax  = (
	    Field((p.n[1],), ((0,), (1,))),
	    Field((p.n[2],), ((0,), (1,)))
	)
	
	assign!(ax[1], fieldgen(i -> i*p.h - 0.5))
	assign!(ax[2], fieldgen(i -> i*p.h - 0.5))
	
	test_poisson(p, ax)
	# test_elastic(p, ax)
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end