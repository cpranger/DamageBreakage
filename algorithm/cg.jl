using OffsetArrays

function lanczos!(A; h, maxit, minit = 10, quiet = false, λtol)
	(v_0, v_1, v_2, u) = h

	γ = OffsetArray(zeros(maxit),    1:maxit)
	η = OffsetArray(zeros(maxit+1),  1:maxit+1)
	λ = OffsetArray(zeros(2),  -1:0)
	Λ = OffsetArray(zeros(2),  -1:0)
	
	assign!(v_0, 0)
	if A(v_0) isa Tuple && !(v_0 isa Tuple) # we have a BC
		assign!(v_1, (fieldgen((_...) -> rand()), A(v_0)[2]))
	else
		assign!(v_1, fieldgen((_...) -> rand()))
	end
	assign!(v_1, v_1 / l2(v_1))
    
	η[1] = 0
	
	# From Meurant et al (2.1)
	k = 1
	for outer k in 1:maxit
		assign!(u, A(v_1) - η[k] * v_0)
		γ[k] = dot(u, v_1)
		assign!(v_2, u - γ[k] * v_1)
		η[k+1] = l2(v_2); η[k+1] > 10*eps(η[k+1]) || break
		assign!(v_2, v_2 / η[k+1])
		assign!(v_0, v_1)
		assign!(v_1, v_2)
		
		# println("(γ, η) = ($(γ[k]), $(η[k+1]))")

		T = LinearAlgebra.SymTridiagonal(γ[1:k], η[2:k])
		(Λ[-1], λ[-1]) = (Λ[0], λ[0])
		(Λ[ 0], λ[ 0]) = extrema <| LinearAlgebra.eigvals(T)
		(λres, Λres)   = map(l -> abs(l[0] - l[-1])/abs(l[0]), (λ, Λ))
		!quiet && println("Lanczos: k = $k, λ = $(λ[0]), Λ = $(Λ[0]), λres = $Λres")
		k >= minit && Λres < λtol #= && λres < λtol=# && break
	end
	return (Λ[0], λ[0])
end



pc_jacobi_apply(a::Tuple,  b) = a*b
pc_jacobi_apply(a::Tensor, b) = StaggeredKernels.TensorOp(:*, a, b)
pc_jacobi_apply(a::Field,  b) = StaggeredKernels.ScalarOp(:*, a, b)

function cg_pc_jacobi!(A::AA, x::X, b::B; h::H, rtol, maxit, minit, quiet = false, λtol = rtol/10) where {AA, X, B, H}
	# after http://www.math.iit.edu/~fass/477577_Chapter_16.pdf
	# with embedded Lanczos process from TODO

	(j, r, z, p, Ap) = h
	
	assign!(p, 0*p) # temporary loan of p
	assign!(j, abs <| 1/diag(p, A(p)))
	(M, m) = 1 ./ minmax(j)
	Meta.@show (M, m)
	
	M⁻¹ = x -> pc_jacobi_apply(x, j)
	
	α  = OffsetArray(zeros(2), -2:-1)
	β  = OffsetArray(zeros(2),  -1:0)
	rz = OffsetArray(zeros(2),  -1:0)
	λ  = OffsetArray(zeros(2),  -1:0)
	Λ  = OffsetArray(zeros(2),  -1:0)
	ε  = OffsetArray(zeros(2),   1:2)

	γ  = OffsetArray(zeros(maxit),    1:maxit)
	η  = OffsetArray(zeros(maxit),    2:maxit+1)
	
	α[-1] = 1
	β[ 0] = 0

	assign!(r, b - A(x))
	assign!(z, M⁻¹(r))
	assign!(p, z)
	
	rz[0] = dot(r, z);# rz[0] > atol || return [];
	
	λ[0] = -1
	Λ[0] = -1

	λf  = maxit
	
	k = 1
	for outer k in 1:maxit
		assign!(Ap, A(p))
		pAp = dot(p, Ap)
		abs(pAp) > 10*eps(pAp) || break
		
		# compute step length
		α[-2] =  α[-1]
		α[-1] = rz[ 0] / pAp
		# update the approximate solution
		assign!(x, x + α[-1] * p)
		# update the residual
		assign!(r, r - α[-1] * Ap)
		assign!(z, M⁻¹(r))
		rz[-1] = rz[0]
		rz[ 0] = dot(r, z)
		abs(rz[0]) > 10*eps(rz[0]) || break
		# compute a gradient correction factor
		β[-1] =  β[0]
		β[ 0] = rz[0] / rz[-1]
		# set the new search direction
		assign!(p, z + β[0]*p)
		
		# update tridiagonalization
		γ[k]   = 1 / α[-1] + β[-1] / α[-2]
		η[k+1] = sqrt(β[0]) / α[-1]
		
		η[k+1] == 0 && break
		
		if k <= λf
			T = LinearAlgebra.SymTridiagonal(γ[1:k], η[2:k])
			(Λ[-1], λ[-1]) = (Λ[0], λ[0])
			(Λ[ 0], λ[ 0]) = extrema <| LinearAlgebra.eigvals(T)
			(λres, Λres)   = map(l -> abs(l[0] - l[-1])/abs(l[0]), (λ, Λ))
			if k >= 10 && Λres < λtol#= && λres < λtol=# 
				λf = k
			end
		end
		
		# r = L^T A L L^-1 x - L^T b
		# ε = x - x^*
		# r = A x - L^T b - [A x^* - b = 0]
		#   = A ε
		# z = M^-1 A ε
		# r^T z = (A ε)^T (M^-1 A ε)
		#       = ε^T A^T M^-1 A ε
		#       = ε^T (A^T L) (L^T A) ε
		# sqrt(r^T z) = |AL ε|
		# λ(AL) |ε| <= [sqrt(r^T z) = |AL ε|] <= Λ(AL) |ε|
		# |ε| <= sqrt(r^T z) / λ(AL)
		# σ(AL) in union(σ(A), σ(L))  >>coloquially speaking<<
		# λ(AL) >= λ(A)λ(L)
		# Λ(AL) <= Λ(A)Λ(L)
		# 1/max(λ(A)Λ(L), Λ(A)λ(L)) <= 1/λ(AL) <= 1/λ(A)λ(L)
		# |ε| <= sqrt(r^T z) / λ(AL) <= sqrt(r^T z) / λ(A)λ(L)

		ε[end] = sqrt(rz[0]) / abs(λ[0]*m) / sqrt(m) / l2(x)
		k == 1 && (ε[1] = ε[end])
		
		if  k <= λf || mod(k, 10) == 0
			k <= λf && !quiet && println("cg: k = $k, log10(εr) = $(log10(ε[end])), λ = $(λ[0]*m), Λ = $(Λ[0]*M), λres = $Λres")
			k >  λf && !quiet && println("cg: k = $k, log10(εr) = $(log10(ε[end]))")
		end
		
		k > minit && ε[end] < rtol &&
			(!quiet && println("cg: k = $k, log10(εr) = $(log10(ε[end]))"); break)
	end

	k == maxit && ε[end] >= rtol && !quiet && println("cg failed to converge in $k iterations")

	return (λ[0]*m, Λ[0]*M, ε[end])
end