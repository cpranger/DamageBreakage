#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

s(e, p) = (p.λ_0 * MajorIdentity() + p.μ_0 * MinorIdentity()) * e

F(v, e, p) = (
    v =   divergence(s(e, p)) / p.r_0,
    e =   symgrad(v),
)

function test_hyperbolic(x, y, p, bc)
	v = Vector(p.n, motion_stags, bc.v)
	e = Tensor(p.n, Symmetric, strain_stags)
	assign!(v, sMode(g)[3,3].gen, (p.o, p.n))
	
	assign!(φ, -divergence(grad(f)), (p.o, p.n))
	
	λ_0 = Mode(g)[p.j...].val
	λ = dot(f, φ, (p.o, p.n)) / dot(f, f, (p.o, p.n)) |> abs |> sqrt
	
	println(" λ = $λ")
	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
	assign!(r, λ^2*f + divergence(grad(f)), (p.o, p.n))
	
	println("||λ^2*f + div(grad(f))|| = $(sqrt <| dot(r, r, (p.o, p.n)))")
	
	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "λ^2 f + div(grad(f))", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
end

# function test_parabolic(x, y, p, bc)
# 	v   = Vector(p.n, motion_stags)
# 	w   = Vector(p.n, motion_stags, bc.v)
# 	r1  =  Field(p.n, div_stags,    Natural())
# 	φ   = Vector(p.n, motion_stags, bc.v)
# 	r2  = Vector(p.n, motion_stags, bc.v)
# 	r3  = Vector(p.n, motion_stags, bc.v)
	
# 	assign!(v, sMode(w)[p.j...].gen, (p.o, p.n))
# 	assign!(φ, curl(curl(v)), (p.o, p.n))

# 	λ_0 = sMode(w)[p.j...].val
# 	λ = dot(v, φ, (p.o, p.n)) / dot(v, v, (p.o, p.n)) |> abs |> sqrt
	
# 	println(" λ = $λ")
# 	println("|λ - λ_0|/|λ_0| = $(abs(λ - λ_0)/abs(λ_0))")
	
# 	assign!(r1, divergence(v),  (p.o, p.n))
# 	assign!(r2, λ^2*v - curl(curl(v)),  (p.o, p.n))
# 	assign!(r3, λ^2*v + divergence(grad(v)),  (p.o, p.n))
	
# 	println("||div(v)|| = $(sqrt <| dot(r1, r1, (p.o, p.n)))")
# 	println("||λ^2*v - curl(curl(v))|| = $(sqrt <| dot(r2, r2, (p.o, p.n)))")
# 	println("||λ^2*v +  div(grad(v))|| = $(sqrt <| dot(r3, r3, (p.o, p.n)))")
	
# 	plt1 = heatmap(x, y, v,  "v", c = :davos)
# 	plt2 = heatmap(x, y, r2, "λ^2 v - curl(curl(v))", c = :davos)
# 	plt  = plot(plt1, plt2; layout = (1, 2))
	
# 	display(plt)
	
# 	readline()
# end

function parameters(; nb, j, μ_0, r_0, h_t)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	o    =  n .- n .+ 1              # logical origin
	l    =  (n .- 2) .* h            # physical domain size
	
	λ_0  =  μ_0                      # assuming Poisson's ratio ν = 1/4
	k_0  =  λ_0 + (2/3)*μ_0          # Bulk modulus
	
	c_s  =  sqrt(μ_0 / r_0)          # solenoidal  wave speed (m/s)
	c_p  =  sqrt((λ_0+2*μ_0)/r_0)    # compressive wave speed (m/s)
	
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
		"--j"
			help = "comma-separated list of the mode number in each dimension (-1 for largest mode)"
			required = true
			arg_type = (NTuple{N, Int} where N)
    end
	
	p = parameters(; μ_0 = 1, r_0 = 1, h_t = @zeros(1), parse_args(s)...); pprintln(p)
	
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i), (p.o[1], p.n[1]))
	assign!(y, fieldgen(i -> i), (p.o[2], p.n[2]))
	
	test_hyperbolic(x, y, p, ImpermeableFreeSlip())
	test__parabolic(x, y, p, ImpermeableFreeSlip())
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end