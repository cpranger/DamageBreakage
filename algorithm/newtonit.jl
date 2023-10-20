
function newton_update!(x, dx, f; h = eps(Float32))
	assign!(x, x + dx)
end

function newton_update!(x, dx, s::SchurComplement; h = eps(Float32))
	D = linearize(s.f_y, (x,),  s.y , (), h = h)

	# (A, B; C, D) (dx, dy) = -(f_x(x, y), f_y(x, y))
	# dy  =  D^-1 (-f_y(x, y) - C dx)
	# dy ?≈ -D^-1   f_y(x + dx, y) [= -D^-1 f_y(x, y) if first x <- x + dx]

	assign!(  x,    x + dx)
	assign!(s.y, s.y - 1/D(1/s.f_y(x, s.y)))
end

function newtonit!(f, u, v, r, h; maxit, atol)
	assign!(r, f(u))
	norm = sqrt <| dot(r, r)
	println("Newton i = 0, ||r|| = $norm")
	
    for i in 1:maxit
		A = linearize(f, u)
		# assign!(v, 0)
		# assign!(v, (0, (A(v) + r)[2]))
		(λ, Λ, ε) = cg_pc_jacobi!(A, v, -r; h = h, rtol = 1e-1, minit = 20, maxit = 100)
		
		newton_update!(u, v, f)

		assign!(r, f(u))
		norm = sqrt <| dot(r, r)
		println("Newton i = $i, ||r|| = $norm, λ = ($(λ[end]), $(Λ[end]))")
		
		norm > atol || break
	end
end
