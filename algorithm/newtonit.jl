function newtonit!(F, u, b, (v, bc_expl), (r, bc_impl), ((λ_1, m_1), (λ_n, m_n)), (h, f); bounds, maxit, atol)
	# axes
	x  = Field((bounds[2][1],), ((0,), (1,)))
	y  = Field((bounds[2][2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (1, bounds[2][1]))
	assign!(y, fieldgen(i -> i), (1, bounds[2][1]))
	
	for i in 1:maxit
		chebyshev!(linearize(F(b), u), v, -F(b)(u), (r, bc_impl); λ = (λ_1, λ_n), v = h, bounds = bounds, atol = atol, maxit = 10*max(bounds[2]...))
		
		assign!(u, u + v, bounds)

		λ_n  = powerit!(linearize(F(b), u), 0,   (m_n, bc_expl), (h, bc_expl); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
		λ_1  = powerit!(linearize(F(b), u), λ_n, (m_1, bc_expl), (h, bc_expl); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
		
		λ_n = rayleighquotientit!(linearize(F(b), u), λ_n + λ_1, (m_n, bc_expl), (r, bc_impl), (h, f), (-λ_1, -λ_n); bounds = bounds, atol = 1e-6, maxit = 2, chebymaxit = max(bounds[2]...))
		λ_1 = rayleighquotientit!(linearize(F(b), u), 0,         (m_1, bc_expl), (r, bc_impl), (h, f), (+λ_1, +λ_n); bounds = bounds, atol = 1e-6, maxit = 2, chebymaxit = max(bounds[2]...))
		
		Meta.@show (λ_1, λ_n)
		
		assign!((r, bc_impl), F(b)(u), bounds)
		norm = sqrt <| dot(r, r, bounds)
		println("Newton i = $i, ||r|| = $norm, λ = ($λ_1, $λ_n)")
		
		plt11 = heatmap(x, y, u, "u", c = :davos)
		plt12 = heatmap(x, y, r, "r", c = :davos)
		plt21 = heatmap(x, y, m_1, "m_1", c = :davos)
		plt22 = heatmap(x, y, m_n, "m_n", c = :davos)
		plt   = plot(plt11, plt12, plt21, plt22; layout = (2, 2))
		display(plt)
		
		norm > atol || break
	end

	return (λ_1, λ_n)
end