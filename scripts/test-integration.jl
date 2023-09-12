#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

s(e, p) = (p.λ_0 * MajorIdentity() + p.μ_0 * MinorIdentity()) * e

F(v, e, p) = (
    v =   divergence(s(e, p)) / p.r_0,
    e =   symgrad(v),
)

function test_hyperbolic(x, y, p)
	readline()
end

function test__parabolic(x, y, p)
	readline()
end

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
    end
	
	p = parameters(; μ_0 = 1, r_0 = 1, h_t = @zeros(1), parse_args(s)...); pprintln(p)
	
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i))
	assign!(y, fieldgen(i -> i))
	
	test_hyperbolic(x, y, p)
	test__parabolic(x, y, p)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end