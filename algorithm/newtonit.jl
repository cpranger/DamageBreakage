
function newton_update!(x, dx, f, bounds; h = eps(Float32))
	assign!(x, x + dx, bounds)
end

function newton_update!(x, dx, s::SchurComplement, bounds; h = eps(Float32))
	D = linearize(s.f_y, (x,),  s.y , (), h = h)

	# (A, B; C, D) (dx, dy) = -(f_x(x, y), f_y(x, y))
	# dy  =  D^-1 (-f_y(x, y) - C dx)
	# dy ?≈ -D^-1   f_y(x + dx, y) [= -D^-1 f_y(x, y) if first x <- x + dx]

	assign!(  x,    x + dx, bounds)
	assign!(s.y, s.y - 1/D(1/s.f_y(x, s.y)), bounds)
end

function newtonit!(f, u, v, r, (h_1, h_2, h_3); bounds, maxit, atol)
	assign!(r, f(u), bounds)
	norm = sqrt <| dot(r, r, bounds)
	println("Newton i = 0, ||r|| = $norm")
	
    for i in 1:maxit
		A = linearize(f, u)
		# assign!(v, 0, bounds)
		assign!(v, (0, (A(v) + r)[2]), bounds)
		(λ, Λ, ε) = cg!(A, v, -r; r = h_1, p = h_2, q = h_3, bounds = bounds, rtol = 1e-1, λtol = 1e-1, minit = 100, maxit = 1000)
		
		newton_update!(u, v, f, bounds)

		assign!(r, f(u), bounds)
		norm = sqrt <| dot(r, r, bounds)
		println("Newton i = $i, ||r|| = $norm, λ = ($(λ[end]), $(Λ[end]))")
		
		norm > atol || break
	end
end
