
# export linearize, SchurComplement

struct SchurComplement # eliminates f(x, y)[2] when linearized
	f
	y
end

(s::SchurComplement)(x) = s.f((x, s.y),)[1]

function linearize(s::SchurComplement, x; h = eps(Float32))
	f_x = (x, y) -> s.f((x, y),)[1]
	f_y = (x, y) -> s.f((x, y),)[2]
	
	A = linearize(f_x, (    ),    x , (s.y,), h = h)
	B = linearize(f_x, (  x,),  s.y , (    ), h = h)
	C = linearize(f_y, (    ),    x , (s.y,), h = h)
	D = linearize(f_y, (  x,),  s.y , (    ), h = h)

	return x -> A(x) - B(1/D(1/C(x)))
end

linearize(f, x; h = eps(Float32)) = v -> imag(f(x + h * im * v)) / h
linearize(f, a::Tuple, x, b::Tuple; h = eps(Float32)) = v -> imag(f(a..., x + h * im * v, b...)) / h

function newton_update!(x, dx, f; h = eps(Float32))
	assign!(x, x + dx)
end

function newton_update!(x, dx, s::SchurComplement; h = eps(Float32))
	f_y = (x, y) -> s.f((x, y),)[2]

	D = linearize(f_y, (x,),  s.y , (), h = h)

	# (A, B; C, D) (dx, dy) = -(f_x(x, y), f_y(x, y))
	# dy  =  D^-1 (-f_y(x, y) - C dx)
	# dy ?≈ -D^-1   f_y(x + dx, y) [= -D^-1 f_y(x, y) if first x <- x + dx]

	dy = -1/D(1/f_y(x, s.y))

	assign!(  x,   x + dx)
	assign!(s.y, s.y + dy)
end

function newtonit!(f, u, r, h; maxit, rtol)
	v = h[1]
	newton_update!(u, 0*u, f)
	
	assign!(r, f(u))
	norm0 = norm = sqrt <| dot(r, r)
	println("Newton i = 0, ||r|| = $norm")
	
    for i in 1:maxit
		A = linearize(f, u)
		
		(_, _, ε) = cg_pc_jacobi!(A, v, -r; h = h[2:end], rtol = 1e-1, minit = 20, maxit = 100)
		
		newton_update!(u, v, f)

		assign!(r, f(u))
		rnorm = l2(r) / norm0
		xnorm = l2(v) / l2(u)
		println("Newton log10‖r_$(i)‖/‖r_0‖ = $(log10(rnorm)), log10‖Δx_$(i)‖/‖x_$(i+1)‖ = $(log10(xnorm))")
		
		rnorm < rtol && xnorm < rtol && break
	end
end
