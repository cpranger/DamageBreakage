cheby(i, x::Complex) = real(0.5*((x + sqrt(x^2 - 1))^i + (x - sqrt(x^2 - 1))^i))

function chebyshev!(A, x, b, (r, bc); λ = (λ_1, λ_2), v, atol, maxit, quiet = false)
	# based on Gutknecht & Röllin (2002; Parallel Computing), algorithm 5.
	assign!((r, bc), b - A(x))
	
	# λ_1 = α - c
	# λ_n = α + c
	α = (λ[end] + λ[1]) / 2
	c = (λ[end] - λ[1]) / 2
	
	ε1 = sqrt <| dot(r, r)
	ε1 > atol || return [ε1]
	
	if abs(c/α) < atol/ε1
		ρ = _ -> abs(c/α)
	else
		ρ = i -> 1/cheby(i,α/c + 0im)
	end

	# Meta.@show α, c, α/c
	# display(plot([ρ(i) for i in 1:100]))
	# readline()

	# massive overstimate, as if simple Jacobi iteration
	nit = min(maxit, max(2,ceil(UInt64, log(atol/ε1)/log(ρ(1)))))
	
	εr = zeros(nit+1)
	ε  = zeros(nit+1)
	ε[1] = ε1
	εr[1] = ε1

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
		
		assign!(v, r - ψ*v)
		assign!(x, x + ω*v)
		assign!((r, bc), r - ω*A(v))
		
		εr[1+i] = sqrt <| dot(r, r)
		ε[1+i]  = ε[1]*ρ(i)
		ε[1+i]  > atol || break

		mod(i, 10) == 1 && !quiet &&
			println("chebyshev: i = $i, log10(ε) = $(log10(εr[1+i])), log10(ε1*ρ) = $(log10(ε[1]*ρ(i)))")
		i += 1
	end

	i == nit+1 && !quiet && println("chebyshev failed to converge in $nit iterations")

	resize!(εr, i); return εr
end