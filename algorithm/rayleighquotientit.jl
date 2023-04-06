rayleigh_quotient(f1, f2, bounds) = dot(f1, f2, bounds) / dot(f1, f1, bounds)

function rayleighquotientit!(A, (m_1, m_n, bc_m), (r, bc_r), (h1, h2), (λ_1, λ_n); bounds, atol, maxit)
	C(λ) = x -> A(x) - λ*x
	
	for i in 1:maxit
		assign!(h1, m_1, bounds)
		ε_1 = chebyshev!(C(0), h1, m_1, (r, bc_r); λ = (λ_1, λ_n), v = h2, bounds = bounds, atol = 0.1, maxit = max(bounds[2]...))
		assign!(m_1, h1 / sqrt(dot(h1, h1, bounds)), bounds)
		assign!((r, bc_r), A(m_1), bounds)
		λ_1 = rayleigh_quotient(m_1, r, bounds)
		assign!(r, r/λ_1 - m_1, bounds)
		r_1 = sqrt <| dot(r, r, bounds)
		
		assign!(h1, m_n, bounds)
		ε_n = chebyshev!(C(λ_n + λ_1), h1, m_n, (r, bc_r); λ = (-λ_1, -λ_n), v = h2, bounds = bounds, atol = 0.1, maxit = max(bounds[2]...))
		assign!(m_n, h1 / sqrt(dot(h1, h1, bounds)), bounds)
		assign!((r, bc_r), A(m_n), bounds)
		λ_n = rayleigh_quotient(m_n, r, bounds)
		assign!(r, r/λ_n - m_n, bounds)
		r_n = sqrt <| dot(r, r, bounds)
		
		display(plot([log10.(ε_1), log10.(ε_n)]))
		# readline()
		
		println("rayleigh $i: λ = ($λ_1, $λ_n), r = ($r_1, $r_n)")

		r_1 > atol || r_n > atol || break
	end
	return (λ_1, λ_n)
end