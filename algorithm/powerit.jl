rayleigh_quotient(f1, f2, bounds) = dot(f1, f2, bounds) / dot(f1, f1, bounds)

function powerit!(op, f1, f2; bounds, atol, maxit)
	fib0 = 1
	fib1 = 1
	fib2 = 2
	iter = 1
    
    norm = sqrt(dot(f1, f1, bounds))
    assign!(f1, f1 / norm, bounds)
    λ = 1.
	
	for i in 1:maxit
		assign!(f2, op(f1)/λ, bounds)
		assign!(f1, op(f2)/λ, bounds)

		if i == fib1
			λ0 = λ; λ = λ*rayleigh_quotient(f2, f1, bounds)

			assign!(f2, (λ0/λ)*f1 - f2, bounds) # -> op(f2)/λ - f2

			resd = sqrt(dot(f2, f2, bounds))
			norm = sqrt(dot(f1, f1, bounds))

			assign!(f1, f1 / norm, bounds) # -> eigenmode
			
			println("powerit: I = $iter, i = $i, Δi = $(fib1-fib0), λ = $λ, |r| = $resd, |f| = $norm")
			
			resd > atol || return λ
			
			fib3 = fib1 + fib2
			fib0 = fib1
			fib1 = fib2
			fib2 = fib3

			iter += 1
		end
	end
	println("power iteration failed to converge in $maxit iterations")
	return λ
end