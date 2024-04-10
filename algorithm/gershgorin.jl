# assumes block matrix in which each block is diagonal
# used for computing the eigenvalues of the explicit submatrix
function gershgorin!(A, h) # three vectors needed
	(x, o, r) = h
	
	assign!(x, 0*x)
	
	result = []

	indices = eachindex(x)
	for self in indices
		others = filter(p -> p != self, indices)

		assign!(x[self], 1)
		assign!(o[self], A(x)[self])
		assign!(x[self], 0)
		
		for other in others
			assign!(x[other], 1)
			assign!(r[self], abs(A(x)[self]))
			assign!(x[other], 0)
		end
		l = min(o[self] - r[self])
		u = max(o[self] + r[self])

		push!(result, (; o = (u + l)/2, r = (u - l)/2))
	end
	return result
end