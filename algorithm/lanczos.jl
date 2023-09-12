import LinearAlgebra
using OffsetArrays

function lanczos!(A, v_1; v_0, v_2, u, maxit)
	α = OffsetArray(zeros(maxit),    1:maxit)
	η = OffsetArray(zeros(maxit+1),  1:maxit+1)

	assign!(v_0, 0)
	assign!(v_1, (fieldgen((_...) -> rand()), A(v_0)[2]))
	assign!(v_1, v_1 / sqrt(dot(v_1, v_1)))
    
	η[1] = 0
	
	# From Meurant et al (2.1)
	
	k = 1
	for outer k in 1:maxit
		assign!(u, A(v_1) - η[k]*v_0)
		α[k] = dot(u, v_1)
		assign!(v_2, u - α[k]*v_1)
		η[k+1] = sqrt <| dot(v_2, v_2);
		η[k+1] > 10*eps(η[k+1]) || break
		assign!(v_2, v_2/η[k+1])
		assign!(v_0, v_1)
		assign!(v_1, v_2)
		
		# println("(α, η) = ($(α[k]), $(η[k+1]))")

		T = LinearAlgebra.SymTridiagonal(α[1:k], η[2:k])
		# Meta.@show T
		λ = extrema <| LinearAlgebra.eigvals(T)
		println("Lanczos $k: $λ")
	end

	T = LinearAlgebra.SymTridiagonal(α[1:k], η[2:k])
	return extrema <| LinearAlgebra.eigvals(T) # Assumes negative-(semi)definite A
end