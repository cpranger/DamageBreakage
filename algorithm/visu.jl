function heatmap(x::Field{Sx}, y::Field{Sy}, f::Field{Sf}, name::String; args...) where {Sx, Sy, Sf}
	stag_f = maximum(Sf) # pick out most central stag
	
	stag_x = (stag_f[1],)
	stag_y = (stag_f[2],)
	
	stag_x in Sx || error("Stag $stag_x not found in x-axis")
	stag_y in Sy || error("Stag $stag_y not found in y-axis")
	
	bounds_f = map((n, s) -> 1:(n - 1 + mod(s, 2)), size(f.data)[2:end], stag_f)
	bounds_x = (bounds_f[1],)
	bounds_y = (bounds_f[2],)
	
	data_x = view(x.data, StaggeredKernels.stagindex(Sx, stag_x), bounds_x...)
	data_y = view(y.data, StaggeredKernels.stagindex(Sy, stag_y), bounds_y...)
	data_f = view(f.data, StaggeredKernels.stagindex(Sf, stag_f), bounds_f...)
	
	return Plots.heatmap(data_x, data_y, data_f'; 
		xlims=(data_x[1], data_x[end]),
		ylims=(data_y[1], data_y[end]),
		# xlabel = "x",
		# ylabel = "y",
		xtickfont=font(3), 
		ytickfont=font(3), 
		guidefont=font(3), 
		legendfont=font(3),
		aspectratio = :equal,#(data_y[end] - data_y[1]) / (data_x[end] - data_x[1]),
		framestyle  = :box,
		# title = name,
		# show = true,
		args...
	)
end

function heatmap(x, y, t::Tensor, name::String; args...)
	plots = map((n, f) -> heatmap(x, y, f, "$name.$n"; args...), keys(t.cpnts), values(t.cpnts))
	return plot(plots...; layout = (length(plots), 1))
end