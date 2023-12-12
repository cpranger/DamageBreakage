# assumes block matrix in which each block is diagonal
# used for computing the eigenvalues of the explicit submatrix
function gershgorin!(A, x::Tuple; h::Tuple{Tuple}) # one help vector needed, supply as tuple for consistency
	d = h[1]
	
	assign!(x, 0*x)
	
	Λ = Inf64
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
		Λ = min(Λ, min(d))
	end
	return Λ
end