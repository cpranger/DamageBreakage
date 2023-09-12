#!/usr/bin/env -S julia --project

include("./header.jl")

using StaggeredKernels.Plane

function test_curl_grad(x, y, p, bc)
	f =  Field(p.n,  div_stags, bc)
	r = Vector(p.n, curl_stags, (z = Natural(),))
	
	assign!(f, Mode(f)[p.j...].gen)
	
	assign!(r, curl(grad(f)))
	
	println("‖curl(grad(f))‖ = ", dot(r, r) |> sqrt)
	
	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
    plt  = plot(plt1, plt2; layout = (2, 1))
	
	display(plt)
	
	readline()
end

function test_div_curl(x, y, p, bc)
	f = Vector(p.n, curl_stags, (z = bc,))
	r =  Field(p.n,  div_stags, Natural())
	
	assign!(f.z, Mode(f.z)[p.j...].gen)
	
	assign!(r, divergence(curl(f)))
	
	println("‖div(curl(f))‖ = ", dot(r, r) |> sqrt)
	
	plt1 = heatmap(x, y, f, "f", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
    plt  = plot(plt1, plt2; layout = (2, 1))
	
	display(plt)
	
	readline()
end

function test_laplacians(x, y, p, bc)
	# div(grad(v)) - grad(div(v)) + curl(curl(v)) = 0
	v = Vector(p.n, motion_stags, bc.v)
	r = Vector(p.n, motion_stags, bc.v)
	
	assign!(v, pMode(v)[p.j...].gen + sMode(v)[p.j...].gen)
	
	assign!(r, grad(divergence(v)) - curl(curl(v)) - divergence(grad(v)))
	
	println("‖grad(div(v)) - curl(curl(v)) - div(grad(v))‖ = ", dot(r, r) |> sqrt)
	
	plt1 = heatmap(x, y, v, "v", c = :davos)
	plt2 = heatmap(x, y, r, "r", c = :davos)
	plt  = plot(plt1, plt2; layout = (1, 2))
	
	display(plt)
	
	readline()
end

function parameters(; nb, j, l_x, h_t)
	n    =  nb .* BLOCK_SIZE         # mesh resolution
	o    =  n .- n .+ 1              # logical origin
	h    =  l_x / (n[1] - 2)         # cell size (square)
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
		"--j"
			help = "comma-separated list of the mode number in each dimension (-1 for largest mode)"
			required = true
			arg_type = (NTuple{N, Int} where N)
    end
	
	p = parameters(; l_x = 1, h_t = @zeros(1), parse_args(s)...); pprintln(p)
	
	# axes
	x  = Field((p.n[1],), ((0,), (1,)))
	y  = Field((p.n[2],), ((0,), (1,)))
	
	assign!(x, fieldgen(i -> i))
	assign!(y, fieldgen(i -> i))
	
	test_curl_grad(x, y, p, Natural())
	test_curl_grad(x, y, p, Essential())
	
	test_div_curl(x, y, p, Natural())
	test_div_curl(x, y, p, Essential())
	
	test_laplacians(x, y, p, PermeableNoSlip())
	test_laplacians(x, y, p, ImpermeableFreeSlip())
	
	"finished!"
end

# see https://stackoverflow.com/a/63385854
if !isdefined(Base, :active_repl)
    main()
end