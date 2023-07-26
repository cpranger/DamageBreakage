import LinearAlgebra
using OffsetArrays

function lanczos!(A, v_1; v_0, v_2, u, bounds, maxit)
	α = OffsetArray(zeros(maxit),    1:maxit)
	η = OffsetArray(zeros(maxit+1),  1:maxit+1)

	assign!(v_0, 0, bounds)
	assign!(v_1, (fieldgen((_...) -> rand()), A(v_0)[2]), bounds)
	assign!(v_1, v_1 / sqrt(dot(v_1, v_1, bounds)), bounds)
    
	η[1] = 0
	
	# From Meurant et al (2.1)
	
	k = 1
	for outer k in 1:maxit
		assign!(u, A(v_1) - η[k]*v_0, bounds)
		α[k] = dot(u, v_1, bounds)
		assign!(v_2, u - α[k]*v_1, bounds)
		η[k+1] = sqrt <| dot(v_2, v_2, bounds);
		η[k+1] > 10*eps(η[k+1]) || break
		assign!(v_2, v_2/η[k+1], bounds)
		assign!(v_0, v_1, bounds)
		assign!(v_1, v_2, bounds)
		
		# println("(α, η) = ($(α[k]), $(η[k+1]))")

		T = LinearAlgebra.SymTridiagonal(α[1:k], η[2:k])
		# Meta.@show T
		λ = extrema <| LinearAlgebra.eigvals(T)
		println("Lanczos $k: $λ")
	end

	T = LinearAlgebra.SymTridiagonal(α[1:k], η[2:k])
	return extrema <| LinearAlgebra.eigvals(T) # Assumes negative-(semi)definite A
end

# function lanczos!(A, v0, v1, w; bounds, maxit)
# 	norm = sqrt <| dot(v0, v0, bounds)
#     assign!(v0, v0 / norm, bounds)
    
# 	α = zeros(0)
# 	β = zeros(0)

# 	# From https://en.wikipedia.org/wiki/Lanczos_algorithm
	
# 	assign!(w, A(v0), bounds)
# 	push!(α, dot(w, v0, bounds))
# 	assign!(w, w - α[end]*v0, bounds)

# 	for j in 2:maxit
# 		β0 = sqrt <| dot(w, w, bounds)
# 		if β0 < 1e-7; break end
# 		push!(β, β0)
# 		assign!(v1, w / β[end], bounds)
		
# 		assign!(w, A(v1), bounds)
# 		push!(α, dot(w, v1, bounds))
		
# 		assign!(w, w - α[end]*v1 - β[end]*v0, bounds)
# 		assign!(v0, v1, bounds)
		
# 		println("(α, η) = ($(α[end]), $(β[end]))")

# 		T = LinearAlgebra.SymTridiagonal(α, β)
# 		# Meta.@show T
# 		λ = extrema <| LinearAlgebra.eigvals(T)
# 		# println("Lanczos $j: $λ")
# 	end

# 	T = LinearAlgebra.SymTridiagonal(α, β)
# 	return extrema <| LinearAlgebra.eigvals(T) # Assumes negative-(semi)definite A
# end