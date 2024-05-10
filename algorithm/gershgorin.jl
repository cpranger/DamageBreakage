
(gershgorin!(A, h::NTuple{N,T}) where {N, T <: AbstractObject}) = gershgorin!(A, map(h -> (h,), h))

# assumes block matrix in which each block is diagonal
# used for computing the eigenvalues of the explicit submatrix
# TODO: can be made about twice as efficient if reductions are implemented on tensor expressions.
function gershgorin!(A, h::NTuple{N,T}) where {N, T <: Tuple} # three vectors needed
	(x, l, u) = h
	
	assign!(x, 0*x)
	
	result = []

	indices = eachindex(x)
	for self in indices
		others = filter(p -> p != self, indices)

		assign!(x[self], 1)
		assign!(l[self], A(x)[self])
		assign!(u[self], A(x)[self])
		assign!(x[self], 0)
		
		for other in others
			assign!(x[other], 1)
			assign!(l[self], l[self] - abs(A(x)[self]))
			assign!(u[self], u[self] - abs(A(x)[self]))
			assign!(x[other], 0)
		end
		push!(result, (; o = (max(l[self]) + min(l[self]))/2, r = (max(l[self]) - min(l[self]))/2))
	end
	return result
end