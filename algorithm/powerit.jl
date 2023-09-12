rayleigh_quotient(f1, f2) = dot(f1, f2) / dot(f1, f1)

(Base.:-(f1::NamedTuple{K}, f2::NamedTuple{K}) where K) = (; zip(K, getfield(f1, k) - getfield(f2, k) for k in K)...)
(Base.:/(f::NamedTuple{K}, d::Number) where K) = (; zip(K, getfield(f, k) / d for k in K)...)
(Base.:*(d::Number, f::NamedTuple{K}) where K) = (; zip(K, d * getfield(f, k) for k in K)...)

function powerit!(A, λ0, (f1, f2, f_bc); atol, maxit)
	fib0 = 1
	fib1 = 1
	fib2 = 2
	iter = 1
    
    B = x -> A(x) - λ0*x
	
	norm = sqrt(dot(f1, f1))
    assign!(f1, f1 / norm)
    λ = 1.
	
	for i in 1:maxit
		assign!((f2, f_bc), B(f1)/λ)
		assign!((f1, f_bc), B(f2)/λ)

		if i == fib1
			λ_ = λ; λ = λ*rayleigh_quotient(f2, f1)

			assign!(f2, (λ_/λ)*f1 - f2) # -> A(f2)/λ - f2

			resd = sqrt(dot(f2, f2))
			norm = sqrt(dot(f1, f1))

			assign!(f1, f1 / norm) # -> eigenmode
			
			println("powerit: I = $iter, i = $i, Δi = $(fib1-fib0), λ = $λ, |r| = $resd, |f| = $norm")
			
			resd > atol || return λ + λ0
			
			fib3 = fib1 + fib2
			fib0 = fib1
			fib1 = fib2
			fib2 = fib3

			iter += 1
		end
	end
	println("power iteration failed to converge in $maxit iterations")
	
	return λ + λ0
end