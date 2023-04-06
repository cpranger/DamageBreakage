import LinearAlgebra

function lanczos!(A, v0, v1, (w, bc); bounds, maxit)
	norm = sqrt <| dot(v0, v0, bounds)
    assign!(v0, v0 / norm, bounds)
    
	α = zeros(0)
	β = zeros(0)

	# From https://en.wikipedia.org/wiki/Lanczos_algorithm
	
	assign!((w, bc), A(v0), bounds)
	push!(α, dot(w, v0, bounds))
	assign!(w, w - α[end]*v0, bounds)

	for j in 2:maxit
		β0 = sqrt <| dot(w, w, bounds)
		if β0 < 1e-7; break end
		push!(β, β0)
		assign!(v1, w / β[end], bounds)
		
		assign!((w, bc), A(v1), bounds)
		push!(α, dot(w, v1, bounds))
		
		assign!(w, w - α[end]*v1 - β[end]*v0, bounds)
		assign!(v0, v1, bounds)
		
		T = LinearAlgebra.SymTridiagonal(α, β)
		λ = extrema <| LinearAlgebra.eigvals(T)
		println("Lanczos $j: $λ")
	end

	T = LinearAlgebra.SymTridiagonal(α, β)
	return extrema <| LinearAlgebra.eigvals(T) # Assumes negative-(semi)definite A
end