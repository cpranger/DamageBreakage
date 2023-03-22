#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

rayleigh_quotient(f1, f2, bounds) = dot(f1, f2, bounds) / dot(f1, f1, bounds)

function powerit!(f1, f2, r, op, bounds, x, y)
	c = Vector(bounds[2], curl_stags, (z = Natural(),))
	
	λ = 1.
	rtol = 1e-6
	niter = 100#00
	fib0 = 1
	fib1 = 1
	fib2 = 2
	iter = 1
	for i in 1:niter
		assign!(f2, op(f1)/λ^2, bounds)
		assign!(f1, op(f2)/λ^2, bounds)

		if i == fib1
			# plt1 = heatmap(x, y, f1, "f1", c = :davos)
			# plt2 = heatmap(x, y, f2, "f2", c = :davos)
			# plt  = plot(plt1, plt2; layout = (1, 2))
			# display(plt)
			
			assign!(r, (f2 + f1)/f2, bounds)

			λ = λ^2*rayleigh_quotient(f2, f1, bounds) |> abs |> sqrt

			rres = sqrt(dot(r, r, bounds))
			norm = sqrt(dot(f1, f1, bounds))
			
			println("powerit: I = $iter, i = $i, Δi = $(fib1-fib0), λ = $λ, |rres| = $rres, |f| = $norm")
			
			if rres < rtol; return λ end

			assign!(f1, f1 / norm, bounds)
			
			fib3 = fib1 + fib2
			fib0 = fib1
			fib1 = fib2
			fib2 = fib3

			iter += 1
		end
	end
	println("power iteration failed to converge in $niter iterations")
	return λ
end

function test_mode(x, y, p, bc)
	f = Field(p.n, div_stags, bc)
	g = Field(p.n, div_stags, bc)
	r = Field(p.n, div_stags, bc)
	
	assign!(f, fieldgen((_...) -> rand()), (p.o, p.n))
	
	λ = 1
	@profview (λ = powerit!(f, g, r, x -> divergence(grad(x)), (p.o, p.n), x, y))
	λ_0 = Mode(g)[-1, -1].val
	
	println(" λ = $λ")
	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
end

function test_s_mode(x, y, p, bc)
	v   = Vector(p.n, motion_stags, bc.v)
	w   = Vector(p.n, motion_stags, bc.v)
	r   = Vector(p.n, motion_stags, bc.v)
	
	assign!(v, (
		x = fieldgen((_...) -> rand()),
		y = fieldgen((_...) -> rand())
	), (p.o, p.n))

	# assign!(w, -curl(curl(v)), (p.o, p.n))
	
	# plt1 = heatmap(x, y, v, "v", c = :davos)
	# plt2 = heatmap(x, y, w, "w", c = :davos)
	# plt  = plot(plt1, plt2; layout = (1, 2))
	# display(plt)
	# readline();
	
	λ = 1
	@profview (λ = powerit!(v, w, r, x -> -curl(curl(x)), (p.o, p.n), x, y))
	λ_0 = sMode(w)[-1, -1].val
	
	println(" λ = $λ")
	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	plt1 = heatmap(x, y, v, "v", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()

	# assign!(v, sMode(w)[p.j...].gen, (p.o, p.n))
	# assign!(φ, curl(curl(v)), (p.o, p.n))

	# λ_0 = sMode(w)[p.j...].val
	# λ = dot(v, φ, (p.o, p.n)) / dot(v, v, (p.o, p.n)) |> abs |> sqrt
	
	# println(" λ = $λ")
	# println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	# assign!(r1, divergence(v),  (p.o, p.n))
	# assign!(r2, λ^2*v - curl(curl(v)),  (p.o, p.n))
	# assign!(r3, λ^2*v + divergence(grad(v)),  (p.o, p.n))
	
	# println("||div(v)|| = $(sqrt <| dot(r1, r1, (p.o, p.n)))")
	# println("||λ^2*v - curl(curl(v))|| = $(sqrt <| dot(r2, r2, (p.o, p.n)))")
	# println("||λ^2*v +  div(grad(v))|| = $(sqrt <| dot(r3, r3, (p.o, p.n)))")
	
	# plt1 = heatmap(x, y, v,  "v", c = :davos)
	# plt2 = heatmap(x, y, r2, "λ^2 v - curl(curl(v))", c = :davos)
	# plt  = plot(plt1, plt2; layout = (1, 2))
	
	# display(plt)
	
	# readline()
end

function test_p_mode(x, y, p, bc)
	v  = Vector(p.n, motion_stags)
	w  = Vector(p.n, motion_stags, bc.v)
	r1 = Vector(p.n, curl_stags, (z = Natural(),))
	φ  = Vector(p.n, motion_stags, bc.v)
	r2  = Vector(p.n, motion_stags, bc.v)
	r3  = Vector(p.n, motion_stags, bc.v)
	
	assign!(v, pMode(w)[p.j...].gen, (p.o, p.n))
	assign!(φ, -grad(divergence(v)), (p.o, p.n))

	λ_0 = pMode(w)[p.j...].val
	λ = dot(v, φ, (p.o, p.n)) / dot(v, v, (p.o, p.n)) |> abs |> sqrt
	
	println(" λ = $λ")
	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	assign!(r1, curl(v), (p.o, p.n))
	assign!(r2, λ^2*v + grad(divergence(v)),  (p.o, p.n))
	assign!(r3, λ^2*v + divergence(grad(v)),  (p.o, p.n))
	
	println("||curl(v)|| = $(sqrt <| dot(r1, r1, (p.o, p.n)))")
	println("||λ^2 v + grad(div(v))|| = $(sqrt <| dot(r2, r2, (p.o, p.n)))")
	println("||λ^2 v + div(grad(v))|| = $(sqrt <| dot(r3, r3, (p.o, p.n)))")
	
	plt1 = heatmap(x, y, v,  "v", c = :davos)
	plt2 = heatmap(x, y, r2, "λ^2 v + grad(div(v))", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
end

function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	o    =  n .- n .+ 1              # logical origin
	
	# collect all variables local to this function:
	vars = Base.@locals
	
	# and return as a named tuple
	return NamedTuple{Tuple(keys(vars))}(values(vars))
end

function main()
	s = ArgParseSettings()
	
	@add_arg_table s begin
    	"--nb"
			help = "comma-separated list of the number of blocks in each dimension"
			required = true
			arg_type = (NTuple{N, Int} where N)
    end
	
	p = parameters(; parse_args(s)...); pprintln(p)
	
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	# test_mode(x, y, p, Essential())
	test_s_mode(x, y, p, ImpermeableFreeSlip())
	# test_p_mode(x, y, p, ImpermeableFreeSlip())

	# Profile.print()
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end