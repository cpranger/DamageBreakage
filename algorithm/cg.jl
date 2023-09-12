using OffsetArrays
import LinearAlgebra

function cg!(A, x, b; r, p, q, rtol, maxit, minit, quiet = false, λtol = rtol)
	γ = OffsetArray(zeros(maxit+1), -1:maxit-1)
	β = OffsetArray(zeros(maxit+1),  0:maxit)
	α = OffsetArray(zeros(maxit),    1:maxit)
	η = OffsetArray(zeros(maxit),    2:maxit+1)
	ρ = OffsetArray(zeros(maxit+1),  0:maxit)
	λ = OffsetArray(zeros(maxit+1),  0:maxit)
	Λ = OffsetArray(zeros(maxit+1),  0:maxit)
	ε = OffsetArray(zeros(maxit),    1:maxit)

	γ[-1] = 1
	β[ 0] = 0

	assign!(r, b - A(x))
	assign!(p, r)
	
	ρ[0] = dot(r, r);# ρ[0] > atol || return [];
	
	λ[0] = -1
	Λ[0] = -1

	λf  = maxit
	χ   = dot(x, x)

	k = 1
	for outer k in 1:maxit
		assign!(q, A(p))
		τ = dot(p, q); abs(τ) > 10*eps(τ) || break;
		γ[k-1] = ρ[k-1] / τ
		assign!(x, x + γ[k-1] * p)
		assign!(r, r - γ[k-1] * q)
		ρ[k] = dot(r, r); abs(ρ[k]) > 10*eps(ρ[k]) || break;
		β[k] = ρ[k] / ρ[k-1]
		assign!(p, r + β[k]*p)
		
		α[k]   = 1 / γ[k-1] + β[k-1] / γ[k-2]
		η[k+1] = sqrt(β[k]) / γ[k-1]
		# println("(α, η) = ($(α[k]), $(η[k+1]))")

		η[k+1] == 0 && break
		
		if k <= 10 || k <= λf
			T = LinearAlgebra.SymTridiagonal(α[1:k], η[2:k])
			(Λ[k], λ[k]) = extrema <| LinearAlgebra.eigvals(T)
			if abs((λ[k] - λ[k-1])/λ[k]) < λtol && abs((Λ[k] - Λ[k-1])/Λ[k]) < λtol
				λf = k
			end
		end

		ε[k] = sqrt(ρ[k]) / abs(λ[min(k, λf)]) / sqrt(dot(x, x))

		if mod(k, 10) == 1
			k <= λf && !quiet && println("cg: k = $k, log10(εr) = $(log10(ε[k])), λ = $(λ[min(k, λf)]), Λ = $(Λ[min(k, λf)])")
			k >  λf && !quiet && println("cg: k = $k, log10(εr) = $(log10(ε[k]))")
		end
		
		k > minit && ε[k] < rtol && break
	end

	k == maxit && ε[k] >= rtol && !quiet && println("cg failed to converge in $k iterations")

	return (λ[1:λf], Λ[1:λf], ε[1:k])
end
