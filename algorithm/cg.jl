using OffsetArrays
import LinearAlgebra

function cg!(A, x, b; r, p, q, bounds, atol, maxit, quiet = false)
	γ = OffsetArray(zeros(maxit+1), -1:maxit-1)
	β = OffsetArray(zeros(maxit+1),  0:maxit)
	α = OffsetArray(zeros(maxit),    1:maxit)
	η = OffsetArray(zeros(maxit),    2:maxit+1)
	ρ = OffsetArray(zeros(maxit+1),  0:maxit)

	γ[-1] = 1
	β[ 0] = 0
	λ = (0, 0)

	assign!(r, b - A(x), bounds)
	assign!(p, r,        bounds)
	
	ρ[0] = dot(r, r, bounds); ρ[0] > atol || return [];
	
	k = 1
	for outer k in 1:maxit
		assign!(q, A(p), bounds)
		τ = dot(p, q, bounds); abs(τ) > 10*eps(τ) || break;
		γ[k-1] = ρ[k-1] / τ
		assign!(x, x + γ[k-1] * p, bounds)
		assign!(r, r - γ[k-1] * q, bounds)
		ρ[k] = dot(r, r, bounds); abs(ρ[k]) > 10*eps(ρ[k]) || break;
		β[k] = ρ[k] / ρ[k-1]
		assign!(p, r + β[k]*p, bounds)
		
		α[k]   = 1 / γ[k-1] + β[k-1] / γ[k-2]
		η[k+1] = sqrt(β[k]) / γ[k-1]
		# println("(α, η) = ($(α[k]), $(η[k+1]))")

		η[k+1] == 0 && break
		
		if mod(k, 10) == 1
			T = LinearAlgebra.SymTridiagonal(α[1:k], η[2:k])
			λ = extrema <| LinearAlgebra.eigvals(T)
			!quiet && println("cg: k = $k, log10(ε) = $(log10(sqrt(ρ[k])/abs(λ[2]))), λ = $λ")
		end
		
		sqrt(ρ[k])/abs(λ[2]) < atol && break
	end

	k == maxit && sqrt(ρ[k])/abs(λ[2]) >= atol && !quiet && println("cg failed to converge in $k iterations")

	return sqrt.(ρ[0:k])./abs(λ[2])
end

# function cg!(A, x, b; r, p, Ap, bounds, atol, maxit, quiet = false)
	
# 	assign!(r, b - A(x), bounds)

# 	assign!(p, r, bounds)
# 	assign!(Ap, A(p), bounds)
	
# 	ρ = zeros(0)
# 	push!(ρ, dot(r, r, bounds))

# 	α = zeros(0)
# 	push!(α, ρ[end] / dot(p, Ap, bounds))

# 	d0 = zeros(0)
# 	d1 = zeros(0)
# 	push!(d0, 1/α[end])

# 	λ = (0, 0)

# 	i = 1; while i <= maxit
# 		assign!(x, x + α[end] *  p, bounds)
# 		assign!(r, r - α[end] * Ap, bounds)

# 		push!(ρ, dot(r, r, bounds))
		
# 		assign!(p, r + ρ[end] / ρ[end-1] * p, bounds)
# 		assign!(Ap, A(p), bounds)

# 		push!(α, ρ[end] / dot(p, Ap, bounds))

# 		if i < 10
# 			push!(d0, 1/α[end] + ρ[end]/ρ[end-1]/α[end-1])
# 			push!(d1, -sqrt(ρ[end]/ρ[end-1])/α[end-1])

# 			T = LinearAlgebra.SymTridiagonal(d0, d1)
# 			λ = extrema <| LinearAlgebra.eigvals(T)
# 			println("Lanczos $i: $λ")
# 		end
		
# 		(mod(i, 10) == 1 || ρ[end] < atol^2) && !quiet &&
# 			println("cg: i = $i, log10(ε) = $(.5*log10(ρ[end]))")
		
# 		ρ[end] < atol^2 && break
		
# 		i += 1
# 	end

# 	i == maxit+1 && !quiet && println("cg failed to converge in $maxit iterations")

# 	resize!(ρ, i); return sqrt.(ρ)
# end