cheby(i, x::Complex) = real(0.5*((x + sqrt(x^2 - 1))^i + (x - sqrt(x^2 - 1))^i))

function chebyshev!(A, x, b; λ = (λ_1, λ_2), v, r, bounds, atol)
	# based on Gutknecht & Röllin (2002; Parallel Computing), algorithm 5.
	assign!(r, b - A(x), bounds)
	
	# λ_1 = α - c
	# λ_n = α + c
	α = (λ[end] + λ[1]) / 2
	c = (λ[end] - λ[1]) / 2
	
	ε1 = sqrt <| dot(r, r, bounds)
	ε1 > atol || return [ε1]
	
	if abs(c/α) < atol/ε1
		ρ = _ -> abs(c/α)
	else
		ρ = i -> 1/cheby(i,α/c + 0im)
	end
	
	# massive overstimate, as if simple Jacobi iteration
	nit = max(2,ceil(UInt64, log(atol/ε1)/log(ρ(1))))
	
	ε = zeros(nit+1)
	ε[1] = ε1

	ψ = 0
	ω = 1/α
	
	i = 1; while i <= nit
		if i == 2
			ψ = -(1/2)*(c/α)^2
			ω = 1/(α - c^2/(2*α))
		elseif i >= 3
			ψ =      - (c/2)^2 * ω^2
			ω = 1/(α - (c/2)^2 * ω)
		end
		
		assign!(v, r - ψ*v , bounds)
		assign!(x, x + ω*v , bounds)
		assign!(r, r - ω*A(v), bounds)
		
		ε[1+i] = sqrt <| dot(r, r, bounds)
		ε[1+i] > atol || break

		mod(i, 10) == 1 &&
			println("chebyshev: i = $i, log10(ε) = $(log10(ε[1+i])), log10(ε1*ρ) = $(log10(ε[1]*ρ(i)))")

		i += 1
	end

	i == nit+1 && println("chebyshev failed to converge in $nit iterations")

	resize!(ε, i); return ε
end