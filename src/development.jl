using ArgCheck
using LinearAlgebra
using Polynomials

function rand_gram_matrix_rank(n, r)

	"""random gram matrix with rank at most r, and with great likelihood r
	"""

	function normalized_random_vector(r)
		v = rand(ComplexF64, r)
		1/norm(v) .* v
	end

	generating_vectors = 0.
	generating_vectors = hcat([normalized_random_vector(r) for i = 1:n]...)

	generating_vectors' * generating_vectors


end


function incomplete_rank(A::Matrix; atol = 1000eps())

	"""computes the permanent of a general matrix n*n A of rank r < n
	following  Two Algorithmic Results for the Traveling Salesman Problem
	Alexander I. Barvinok

	"""

	n = size(A,1)
	r = rank(A, atol = atol)

	@argcheck r < n

	qr_dec = qr(A)
	G = qr_dec.Q
	B = qr_dec.R

end

atol = 1000eps()
A = rand_gram_matrix_rank(4,3)

n = size(A,1)
r = rank(A, atol = atol)

@argcheck r < n

qr_dec = qr(A)
G = qr_dec.Q
B = qr_dec.R

j = 1

G[j,:][1:r]

L = prod([Polynomial(B[:,j][1:r]) for j in 1:n]).coeffs
R = prod([Polynomial(G[i,:][1:r]) for i in 1:n]).coeffs

[Polynomial(B[:,j][1:r]) for j in 1:n]
