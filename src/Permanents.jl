module Permanents

export naive, naive_tensor, ryser

using Combinatorics
using LinearAlgebra
using Test

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

function multi_dim_ryser(U, gram_matrix)

    # https://arxiv.org/pdf/1410.7687.pdf

    n = size(U)[1]
    nstring = collect(1:n)
    sub_nstring = collect(powerset(nstring))
    sub_nstring = sub_nstring[2:length(sub_nstring)]
    res = 0

    function delta_set(S1, S2)
        if S1 == S2
            return 1
        else
            return 0
        end
    end

    for r = 1:length(sub_nstring)
        for s = r:length(sub_nstring)
            R = sub_nstring[r]
            S = sub_nstring[s]

            t = prod(
                sum(
                    U[ss, j] * conj(U[rr, j]) * gram_matrix[rr, ss] for rr in R for ss in S
                ) for j = 1:n
            )
            res +=
                (2 - delta_set(S, R)) * (-1)^(length(S) + length(R)) * real(t)
        end
    end
    return res
end

function ryser(A::AbstractMatrix)
	"""computes the permanent of A using ryser with Gray ordering"""
		# code from https://discourse.julialang.org/t/matrix-permanent/10766
		# see Combinatorial Algorithms for Computers and Calculators, Second Edition (Computer Science and Applied Mathematics) by Albert Nijenhuis, Herbert S. Wilf (z-lib.org)
		# chapter 23 for the permanent algorithm
		# chapter 1 for the gray code

	    function grayBitToFlip(n::Int)
	    	n1 = (n-1) ⊻ ((n-1)>>1)
	    	n2 = n ⊻ (n>>1)
	    	d = n1 ⊻ n2

	    	j = 0
	    	while d>0
	    		d >>= 1
	    		j += 1
	    	end
	    	j
	    end

	    n,m = size(A)
		if (n == m)
	        AA = 2*A
			D = true
			rho = [true for i = 1:n]
			v = sum(A, dims = 2)
			p = prod(v)
			@inbounds for i = 1:(2^(n-1)-1)
				a = grayBitToFlip(i)+1
	            if rho[a]
	                @simd for j=1:n
	                    v[j] -= AA[j,a]
	                end
	            else
	                @simd for j=1:n
	                    v[j] += AA[j,a]
	                end
	            end
	            rho[a] = !rho[a]
	            pv = one(typeof(p))
	            @simd for j=1:n
	                pv *= v[j] #most expensive part
	            end
				D = !D
	            if D
	                p += pv
	            else
	                p -= pv
	            end
			end

			return p * 2.0^(1-n)
		else
			throw(ArgumentError("perm: argument must be a square matrix"))
		end

end

function fast_glynn_perm(U::AbstractMatrix{T}) where T

	""" https://codegolf.stackexchange.com/questions/97060/calculate-the-permanent-as-quickly-as-possible """

	size(U)[1] == size(U)[2] ? n=size(U)[1] : error("Non square matrix as input")

	row_ = [U[:,i] for i in 1:n]
	row_comb = [sum(row) for row in row_]

	res = 0
	old_gray = 0
	sign = +1
	binary_power_dict = Dict(2^k => k for k in 0:n)
	num_iter = 2^(n-1)

	for bin_index in 1:num_iter
		res += sign * reduce(*, row_comb)

		new_gray = bin_index ⊻ trunc(Int, bin_index/2)
		gray_diff = old_gray ⊻ new_gray
		gray_diff_index = binary_power_dict[gray_diff]

		new_vec = U[gray_diff_index+1,:]
		direction = 2 * cmp(old_gray, new_gray)

		for i in 1:n
			row_comb[i] += new_vec[i] * direction
		end

		sign = -sign
		old_gray = new_gray
	end

	return res/num_iter

end

# Function to compute the multinomial coefficient for a given tuple of integers
function multinomial_coeff(alpha::Vector{Int})
    return factorial(sum(alpha)) / prod(factorial(a) for a in alpha)
end

function extract_polynomial_terms(poly)
    # Get variables of the polynomial
    vars = variables(poly)
    
    # Get coefficients and monomials
    coeffs, monos = coefficients(poly), monomials(poly)
    
    # Initialize array to store results
    # Each element will be a tuple (coefficient, exponents)
    terms = []
    
    # Process each term
    for (coeff, mono) in zip(coeffs, monos)
        # Get exponents for this monomial
        exp = exponents(mono)
        
        # Add tuple of (coefficient, exponents) to results
        push!(terms, (coeff, exp))
    end
    
    return terms
end

function polynomial_to_dict(poly)
    # Get coefficients and monomials
    coeffs, monos = coefficients(poly), monomials(poly)
    
    # Create dictionary
    terms_dict = Dict{Vector{Int}, Complex{Float64}}()
    
    # Fill dictionary with terms
    for (coeff, mono) in zip(coeffs, monos)
        exp = exponents(mono)
        terms_dict[exp] = coeff
    end
    
    return terms_dict
end

# Helper function to print the dictionary in a readable format
function print_polynomial_dict(dict)
    println("Terms:")
    for (exp, coeff) in sort(collect(dict), by=first)
        println("$exp => $coeff")
    end
end

Base.factorial(coeff::Vector) = prod(factorial(c) for c in coeff)

function multiply_shared_terms_with_coefficient(dict1::Dict{Vector{Int}, ComplexF64}, 
                             dict2::Dict{Vector{Int}, ComplexF64})
    # Initialize result dictionary
    result = zero(ComplexF64)
    
    # Find shared keys (exponents)
    shared_exponents = intersect(keys(dict1), keys(dict2))
    
    # Multiply coefficients for shared exponents
    for exp in shared_exponents
        
        result += factorial(exp) * dict1[exp] * dict2[exp]
    end
    
    return result
end

# Helper function to print the results
function print_multiplication_results(dict1::Dict{Vector{Int}, ComplexF64}, 
                                   dict2::Dict{Vector{Int}, ComplexF64})
    result = multiply_shared_terms(dict1, dict2)
    
    println("Shared terms and their products:")
    for (exp, coeff) in sort(collect(result), by=first)
        println("\nExponents: $exp")
        println("First coefficient: $(dict1[exp])")
        println("Second coefficient: $(dict2[exp])")
        println("Product: $coeff")
    end
    
    println("\nTotal number of shared terms: $(length(result))")
end


# Function to compute the permanent using Barvinok's approach
function incomplete_rank(A::AbstractMatrix)

	A = Matrix(A)

    n = size(A, 1)
    r = rank(A, atol = atol)
    
    if r == n 
        @warn"while the algorithm will work for any rank, you should rather use another one if full rank, using ryser instead"
        return ryser(A)
    end

    # Perform QR decomposition
    qr_dec = qr(A)
    G = qr_dec.Q
    B = qr_dec.R

    # Define polynomial variables x1, x2, ..., xr
    @polyvar x[1:r]

    # Construct the polynomials L and R using the different variables
    L_poly = prod([sum(B[:, j][i] * x[i] for i in 1:r) for j in 1:n])
    R_poly = prod([sum(G[i, :][j] * x[j] for j in 1:r) for i in 1:n])

    # Extract the coefficients of the polynomials L and R, including multinomial factors
    coeffs_L = polynomial_to_dict(L_poly)
    coeffs_R = polynomial_to_dict(R_poly)

    return multiply_shared_terms_with_coefficient(coeffs_L, coeffs_R)

end

end
