#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

function test_mode(x, y, p, bc)
	f = Field(p.n, div_stags)
	g = Field(p.n, div_stags)
	
	assign!(f, fieldgen((_...) -> rand()), (p.o, p.n))
	
	A = x -> divergence(grad(x))
	λ = powerit!(A, 0, (f, bc), (g, bc); bounds = (p.o, p.n), maxit = 10000, atol = 1e-7)
	λ_0 = -(Mode(f)[-1, -1].val)^2
	
	println(" λ = $λ, λ_0 = $λ_0")
	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, g, "g", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
end

function test_s_mode(x, y, p, bc)
	v   = Vector(p.n, motion_stags)
	w   = Vector(p.n, motion_stags)
	
	assign!(v, (
		x = fieldgen((_...) -> rand()),
		y = fieldgen((_...) -> rand())
	), (p.o, p.n))

	A = x -> -curl(curl(x))
	λ = powerit!(A, 0, (v, bc.v), (w, bc.v); bounds = (p.o, p.n), maxit = 10000, atol = 1e-7)
	λ_0 = -(sMode(v)[-1, -1].val)^2
	
	println(" λ = $λ, λ_0 = $λ_0")
	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	plt1 = heatmap(x, y, v, "v", c = :davos)
	plt2 = heatmap(x, y, w, "w", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
end

function test_p_mode(x, y, p, bc)
	v   = Vector(p.n, motion_stags)
	w   = Vector(p.n, motion_stags)
	
	assign!(v, (
		x = fieldgen((_...) -> rand()),
		y = fieldgen((_...) -> rand())
	), (p.o, p.n))

	A = x -> grad(divergence(x))
	λ = powerit!(A, 0, (v, bc.v), (w, bc.v); bounds = (p.o, p.n), maxit = 10000, atol = 1e-7)
	λ_0 = -(pMode(v)[-1, -1].val)^2
	
	println(" λ = $λ, λ_0 = $λ_0")
	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	plt1 = heatmap(x, y, v, "v", c = :davos)
	plt2 = heatmap(x, y, w, "w", c = :davos)
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
	
	  test_mode(x, y, p, Essential())
	test_s_mode(x, y, p, ImpermeableFreeSlip())
	test_p_mode(x, y, p, ImpermeableFreeSlip())
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end