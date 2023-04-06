#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

using StaggeredKernels.Plane

function test_mode(x, y, p, bc)
	v0 = Field(p.n, div_stags)
	v1 = Field(p.n, div_stags)
	w  = Field(p.n, div_stags)
	
	assign!((v0, bc), fieldgen((_...) -> rand()), (p.o, p.n))
	
	A = x -> divergence(grad(x))
	(λn, λ1) = lanczos!(A, v0, v1, (w, bc); bounds = (p.o, p.n), maxit = max(p.n...))
	λ_0 = -(Mode(w, bc)[-1, -1].val)^2

	println(" λ = $λn, λ_0 = $λ_0")
	println("|λ - λ_0|/|λ_0| = $(abs(λn - λ_0)/abs(λ_0))")
	
	readline()
end

function test_s_mode(x, y, p, bc)
	v0 = Vector(p.n, motion_stags)
	v1 = Vector(p.n, motion_stags)
	w  = Vector(p.n, motion_stags)
	
	assign!((v0, bc.v), (
		x = fieldgen((_...) -> rand()),
		y = fieldgen((_...) -> rand())
	), (p.o, p.n))

	A = x -> -curl(curl(x))
	(λn, λ1) = lanczos!(A, v0, v1, (w, bc.v); bounds = (p.o, p.n), maxit = max(p.n...))
	λ_0 = -(sMode(v0, bc.v)[-1, -1].val)^2
	
	println(" λ = $λn, λ_0 = $λ_0")
	println("|λ - λ_0|/|λ_0| = $(abs(λn - λ_0)/abs(λ_0))")
	
	readline()
end

function test_p_mode(x, y, p, bc)
	v0 = Vector(p.n, motion_stags)
	v1 = Vector(p.n, motion_stags)
	w  = Vector(p.n, motion_stags)
	
	assign!((v0, bc.v), (
		x = fieldgen((_...) -> rand()),
		y = fieldgen((_...) -> rand())
	), (p.o, p.n))

	A = x -> grad(divergence(x))
	(λn, λ1) = lanczos!(A, v0, v1, (w, bc.v); bounds = (p.o, p.n), maxit = max(p.n...))
	λ_0 = -(pMode(v0, bc.v)[-1, -1].val)^2
	
	println(" λ = $λn, λ_0 = $λ_0")
	println("|λ - λ_0|/|λ_0| = $(abs(λn - λ_0)/abs(λ_0))")
	
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