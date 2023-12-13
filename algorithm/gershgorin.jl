# assumes block matrix in which each block is diagonal
# used for computing the eigenvalues of the explicit submatrix
function gershgorin!(A, h) # two vectors needed
	(x, d) = h
	
	assign!(x, 0*x)
	
	Λ = [Inf, -Inf]
	indices = eachindex(x)
	for self in indices
		others = filter(p -> p != self, indices)

		assign!(x[self], 1)
		assign!(d[self], A(x)[self])
		assign!(x[self], 0)
		for other in others
			assign!(x[other], 1)
			assign!(d[self],  d[self] - abs(A(x)[self]))
			assign!(x[other], 0)
		end
		Λ[1] = min(Λ[1], min(min.(d)...))

		assign!(x[self], 1)
		assign!(d[self], A(x)[self])
		assign!(x[self], 0)
		for other in others
			assign!(x[other], 1)
			assign!(d[self],  d[self] + abs(A(x)[self]))
			assign!(x[other], 0)
		end
		Λ[2] = max(Λ[2], max(max.(d)...))
	end
	return Λ
end