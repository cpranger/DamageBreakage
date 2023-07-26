# linearize(F, x; h = eps(Float32)) = v -> imag(F(x + h * im * v)) / h

# linearize(F, a::Tuple, x, b::Tuple; h = eps(Float32)) = v -> imag(F(a..., x + h * im * v, b...)) / h

struct SchurComplement # eliminates CD y = v
	AB
	CD
	y
	v
end

(f::SchurComplement)(x) = f.AB(x, f.y)

# function linearize(f::SchurComplement, x; h = eps(Float32))
# 	A = linearize(f.AB, (    ),    x , (f.y,), h = h)
# 	B = linearize(f.AB, (  x,),  f.y , (    ), h = h)
# 	C = linearize(f.CD, (    ),    x , (f.y,), h = h)
# 	D = linearize(f.CD, (  x,),  f.y , (    ), h = h)

# 	return x -> A(x) - B(1/D(1/C(x)))
# end

function schur_update!(f::SchurComplement, x, dx, bounds; h = eps(Float32))
	C = linearize(f.CD, (    ),    x , (f.y,), h = h)
	D = linearize(f.CD, (  x,),  f.y , (    ), h = h)

	assign!(f.v, -f.CD(x, f.y), bounds)
	assign!(f.y,  f.y + 1/D(1/(f.v - C(x + dx))), bounds)
end

# function schur_update!(f::SchurComplement, x)
# 	assign!(v, -CD(x, f.y))
# 	assign!(f.y,  f.y + 1/D(1/(v - C(x))))
# end

function newtonit!(f, u, (v, bc_expl), (r, bc_impl), ((λ_1, m_1), (λ_n, m_n)), (h1, h2); bounds, maxit, atol)
	# axes
	ax_x  = Field((bounds[2][1],), ((0,), (1,)))
	ax_y  = Field((bounds[2][2],), ((0,), (1,)))
	
	assign!(ax_x, fieldgen(i -> i), (1, bounds[2][1]))
	assign!(ax_y, fieldgen(i -> i), (1, bounds[2][1]))
	
	assign!((r, bc_impl), f(u), bounds)
	norm = sqrt <| dot(r, r, bounds)
	println("Newton i = 0, ||r|| = $norm, λ = ($λ_1, $λ_n)")
	
    for i in 1:maxit
		chebyshev!(linearize(f, u), v, -r, (r, bc_impl); λ = (λ_1, λ_n), v = h1, bounds = bounds, atol = atol, maxit = 10*max(bounds[2]...))
		
		assign!(u, u + v, bounds)
        hasmethod(schur_update!, typeof((f, u, v, bounds))) && schur_update!(f, u, v, bounds)

		λ_n  = powerit!(linearize(f, u), 0,   (m_n, h1, bc_expl); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
		λ_1  = powerit!(linearize(f, u), λ_n, (m_1, h1, bc_expl); bounds = bounds, maxit = max(bounds[2]...), atol = 1e-3)
		
		λ_n = rayleighquotientit!(linearize(f, u), λ_n + λ_1, (m_n, bc_expl), (r, bc_impl), (h1, h2), (-λ_1, -λ_n); bounds = bounds, atol = 1e-6, maxit = 2, chebymaxit = max(bounds[2]...))
		λ_1 = rayleighquotientit!(linearize(f, u), 0,         (m_1, bc_expl), (r, bc_impl), (h1, h2), (+λ_1, +λ_n); bounds = bounds, atol = 1e-6, maxit = 2, chebymaxit = max(bounds[2]...))
		
		Meta.@show (λ_1, λ_n)
		
		assign!((r, bc_impl), f(u), bounds)
		norm = sqrt <| dot(r, r, bounds)
		println("Newton i = $i, ||r|| = $norm, λ = ($λ_1, $λ_n)")
		
		plt11 = heatmap(ax_x, ax_y, u, "u", c = :davos)
		plt12 = heatmap(ax_x, ax_y, r, "r", c = :davos)
		plt21 = heatmap(ax_x, ax_y, m_1, "m_1", c = :davos)
		plt22 = heatmap(ax_x, ax_y, m_n, "m_n", c = :davos)
		plt   = plot(plt11, plt12, plt21, plt22; layout = (2, 2))
		display(plt)
		
		norm > atol || break
	end

	return (λ_1, λ_n)
end
