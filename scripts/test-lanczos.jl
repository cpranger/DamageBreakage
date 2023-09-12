#!/usr/bin/env -S julia --project

using Revise

includet("./header.jl")

gen_ones() = fieldgen((_...) -> 1)
gen_rand() = fieldgen((_...) -> rand())

gen_ones(f::Field) = gen_ones()
gen_rand(f::Field) = gen_rand()

(gen_ones(v::Tensor{S, NamedTuple{N}}) where {S, N}) = (; zip(N, gen_ones() for n in N...)...)
(gen_rand(v::Tensor{S, NamedTuple{N}}) where {S, N}) = (; zip(N, gen_rand() for n in N...)...)

(cstep(f, h)) = imag(f) / h
(cstep(bc::Essential{D}, h) where D) = Essential{D}(cstep(bc.expr, h))
(cstep(  bc::Natural{D}, h) where D) =   Natural{D}(cstep(bc.expr, h))
(cstep(f::Tuple{F, BCs}, h) where {F, BCs <: Tuple}) = (cstep(f[1], h), map(bc -> cstep(bc, h), f[2]))

linearize(f, x; h = eps(Float32)) = v -> cstep(f(x + h * im * v), h)

using StaggeredKernels.Plane

function test_mode(x, y, p)
	v_0 = Field(p.n, div_stags)
	v_1 = Field(p.n, div_stags)
	v_2 = Field(p.n, div_stags)
	u   = Field(p.n, div_stags)
	
	bc = x -> (
		Essential(:-, :y, -x),
		Essential(:-, :x, -x),
		Essential(:+, :x, -x),
		Essential(:+, :y, -x)
	)

	assign!(v_1, gen_rand(v_1))
	# assign!(v_1, #=(=#gen_ones(v_1)#=, bc(0))=#)
	# display <| heatmap(x, y, v_1, "v_1", c = :davos)
	# Meta.@show dot(v_1, v_1)
	# return

	mode = Mode(v_0, bc(0))
	(λn_0, λ1_0) = (-mode[-1,-1].val^2, -mode[1,1].val^2)
	println("λ_* = ($λn_0, $λ1_0)")
	
	A = x -> (divergence(grad(x)), bc(x))
	(λn, λ1) = lanczos!(A, v_1; v_0=v_0, v_2=v_2, u=u, maxit = 100)
	
	println(" λ = ($λn, $λ1)")
	println("|λ - λ_*|/|λ_*| = ($(abs(λn - λn_0)/abs(λn_0)), $(abs(λ1 - λ1_0)/abs(λ1_0)))")
end

function test_s_mode(x, y, p, bc)
	v0 = Vector(p.n, motion_stags)
	v1 = Vector(p.n, motion_stags)
	w  = Vector(p.n, motion_stags)
	
	assign!((v0, bc.v), (
		x = fieldgen((_...) -> rand()),
		y = fieldgen((_...) -> rand())
	))

	A = x -> -curl(curl(x))
	(λn, λ1) = lanczos!(A, v0, v1, (w, bc.v); maxit = max(p.n...))
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
	))

	A = x -> grad(divergence(x))
	(λn, λ1) = lanczos!(A, v0, v1, (w, bc.v); maxit = max(p.n...))
	λ_0 = -(pMode(v0, bc.v)[-1, -1].val)^2
	
	println(" λ = $λn, λ_0 = $λ_0")
	println("|λ - λ_0|/|λ_0| = $(abs(λn - λ_0)/abs(λ_0))")
	
	readline()
end

function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	
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
	
	assign!(x, fieldgen(i -> i))
	assign!(y, fieldgen(i -> i))
	
	test_mode(x, y, p)
	# test_s_mode(x, y, p, ImpermeableFreeSlip())
	# test_p_mode(x, y, p, ImpermeableFreeSlip())
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end