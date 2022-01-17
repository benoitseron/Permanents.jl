using Combinatorics
using LinearAlgebra

const ATOL = 1e4 * eps()

module Permanents

function naive(U::AbstractMatrix)

    """ computes the permanent using the definition as a naive implementation"""

    n = size(U)[1]

    sum_diag(perm, U) = prod(U[i,perm[i]] for i = 1:n)

    return sum(sum_diag(sigma, U) for sigma in permutations(collect(1:n)))

end

function naive_tensor(W::Array)

    if length(size(W)) == 3 && all(size(W) .== size(W)[1])
        n = size(W)[1]

        sum_tensor_diag(sigma,rho, W) = prod(W[i,sigma[i],rho[i]] for i = 1:n)
        return sum(sum_tensor_diag(sigma,rho, W) for sigma in permutations(collect(1:n)) for rho in permutations(collect(1:n)))
    else
        error("tensor permanent implemented only for square 3-indices tensors")
    end

end

function ryser(A::AbstractMatrix)

	"""computes the permanent of A using Ryser's algorithm with Gray ordering"""

	# see Combinatorial Algorithms for Computers and Calculators,
	# Second Edition (Computer Science and Applied Mathematics)
	# by Albert Nijenhuis, Herbert S. Wilf (z-lib.org)
	# chapter 23 for the permanent algorithm
	# chapter 1 for the gray code

    n = length(A[1,:])
    if length(A[1,:]) != length(A[:,1])
        error("matrix not square")
    end

    a = zeros(Int, n)
    f = collect(0:n)
    f[n+1] = n

    j = 0

    total = zero(eltype(A))
    sums = zeros(eltype(A), n)
	one_this_type = one(eltype(A))

    @inbounds for parity = 0 : 2^n-2

        j = f[1]
        f[1] = 0

        f[j+1] = f[j+2]
        f[j+2] = j + 1
        a[j+1] = 1 - a[j+1] # there you find the bit to flip, j+1

        this_product = one_this_type

        if a[j+1] == 1
            for i = 1 : n
                sums[i] +=  A[j+1, i]
				this_product *= sums[i]
            end
        else
            for i = 1 : n
                sums[i] -=  A[j+1, i]
				this_product *= sums[i]
            end
        end

		if parity % 2 == 0
        	total += this_product
		else
			total -= this_product
		end

    end

    (-1)^(n-1) * total

end

end
