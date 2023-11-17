
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

function newtonit!(f, u, v, r, h; maxit, rtol)
	newton_update!(u, 0*u, f)
	
	assign!(r, f(u))
	norm0 = norm = sqrt <| dot(r, r)
	println("Newton i = 0, ||r|| = $norm")
	
    for i in 1:maxit
		A = linearize(f, u)
		
		(_, _, ε) = cg_pc_jacobi!(A, v, -r; h = h, rtol = 1e-1, minit = 20, maxit = 100)
		
		newton_update!(u, v, f)

		assign!(r, f(u))
		rnorm = l2(r) / norm0
		xnorm = l2(v) / l2(u)
		println("Newton ‖r_$(i)‖/‖r_0‖ = $rnorm, ‖Δx_$(i)‖/‖x_$(i+1)‖ = $xnorm")
		
		rnorm < rtol && xnorm < rtol && break
	end
end
