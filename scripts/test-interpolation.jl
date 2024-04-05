#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

function test_interp(p, ax)
	f = Field(p.n, ((1, 1), (0, 0)))
	r = Field(p.n, ((0, 1), (1, 0)))
	
    assign!(f, fieldgen((i,j) -> i*j))
    assign!(r, fieldgen((i,j) -> i*j) - interpolate(f))

	plt1 = heatmap(ax..., f, "f", c = :davos)
	plt2 = heatmap(ax..., r, "r", c = :davos)
    plt  = plot(plt1, plt2; layout = (2, 1))
	display(plt)
end

function test_tensorcontract(p, ax)
	o  =  Field(p.n, div_stags)
	v  = Tensor(p.n, motion_stags)
	e  = Tensor(p.n, Symmetric, strain_stags)
	j1 =  Field(p.n, state_stags)
	j2 =  Field(p.n, state_stags)
	j3 =  Field(p.n, state_stags)
    
    assign!(o, fieldgen((i,j) -> i^3 * j^3))
	assign!(v, grad(o))
	assign!(e, symgrad(v))

	assign!(j1, J1(e))
	assign!(j2, J2(e))
	assign!(j3, J3(e))
	
	plt1 = heatmap(ax..., j1, "j1", c = :davos)
	plt2 = heatmap(ax..., j2, "j2", c = :davos)
	plt3 = heatmap(ax..., j3, "j3", c = :davos)
    plt  = plot(plt1, plt2, plt3; layout = (3, 1))
	display(plt)
end


function parameters(; nb)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	o    =  n .- n .+ 1              # logical origin
	h    =  1 / (n[1] - 2)           # cell size (square)
	l    =  (n .- 2) .* h            # physical domain size
	
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
	
	p = parameters(; parse_args(s)...);# pprintln(p)
	
	# axes
	ax  = (
	    Field((p.n[1],), ((0,), (1,))),
	    Field((p.n[2],), ((0,), (1,)))
	)
	
	assign!(ax[1], fieldgen(i -> i*p.h - 0.5))
	assign!(ax[2], fieldgen(i -> i*p.h - 0.5))
	
	# test_interp(p, ax)
    test_tensorcontract(p, ax)
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end