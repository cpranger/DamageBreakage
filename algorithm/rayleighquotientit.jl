rayleigh_quotient(f1, f2) = dot(f1, f2) / dot(f1, f1)

function rayleighquotientit!(A, λ0, (m, bc_m), (r, bc_r), (h1, h2), (λ_1, λ_n); atol, maxit, chebymaxit, chebyquiet = true)
	C(λ) = x -> A(x) - λ*x
	λ = λ0
	for i in 1:maxit
		assign!(h1, m)
		chebyshev!(C(λ0), h1, m, (r, bc_r); λ = (λ_1, λ_n), v = h2, atol = 0.1, maxit = chebymaxit, quiet = chebyquiet)
		assign!(m, h1 / sqrt(dot(h1, h1)))
		assign!((r, bc_r), A(m))
		λ = rayleigh_quotient(m, r);
		assign!(r, r/λ - m)
		ρ = sqrt <| dot(r, r)
		println("rayleigh $(i): λ = ($λ), r = ($ρ)")
		ρ > atol || break
	end
	return λ
end
