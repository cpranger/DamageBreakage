rayleigh_quotient(f1, f2, bounds) = dot(f1, f2, bounds) / dot(f1, f1, bounds)

function rayleighquotientit!(A, λ0, (m, bc_m), (r, bc_r), (h1, h2), (λ_1, λ_n); bounds, atol, maxit, chebymaxit = max(bounds[2]...), chebyquiet = true)
	C(λ) = x -> A(x) - λ*x
	λ = λ0
	for i in 1:maxit
		assign!(h1, m, bounds)
		chebyshev!(C(λ0), h1, m, (r, bc_r); λ = (λ_1, λ_n), v = h2, bounds = bounds, atol = 0.1, maxit = chebymaxit, quiet = chebyquiet)
		assign!(m, h1 / sqrt(dot(h1, h1, bounds)), bounds)
		assign!((r, bc_r), A(m), bounds)
		λ = rayleigh_quotient(m, r, bounds);
		assign!(r, r/λ - m, bounds)
		ρ = sqrt <| dot(r, r, bounds)
		println("rayleigh $(i): λ = ($λ), r = ($ρ)")
		ρ > atol || break
	end
	return λ
end
