module Permanents

export naive, naive_tensor, ryser, multi_dim_ryser, fast_glynn_perm, glynn, glynn_precision, ryser_tensor

using Combinatorics
using LinearAlgebra
using ArgCheck

function naive(U::AbstractMatrix)

    """ computes the permanent using the definition as a naive implementation"""

    n = size(U)[1]

    sum_diag(perm, U) = prod(U[i,perm[i]] for i = 1:n)

    return sum(sum_diag(sigma, U) for sigma in permutations(collect(1:n)))

end
function is_a_square_three_tensor(W::Array)

	length(size(W)) == 3 && all(size(W) .== size(W)[1])

end


function naive_tensor(W::Array)

    if is_a_square_three_tensor(W)
        n = size(W)[1]

        sum_tensor_diag(sigma,rho, W) = prod(W[i,sigma[i],rho[i]] for i = 1:n)
        return sum(sum_tensor_diag(sigma,rho, W) for sigma in permutations(collect(1:n)) for rho in permutations(collect(1:n)))
    else
        error("tensor permanent implemented only for square 3-indices tensors")
    end

end

function ryser_tensor(W::Array)

	"""implements the tensor permanent according to Ryser's factorization
	without using Gray factoring

	corresponds to Eq. 21 of https://arxiv.org/abs/1410.7687"""

	@argcheck is_a_square_three_tensor(W) "tensor permanent implemented only for square 3-indices tensors"

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
					W[rr,ss,j] for rr in R for ss in S
                ) for j = 1:n
            )
            res +=
                (2 - delta_set(S, R)) * (-1)^(length(S) + length(R)) * real(t)
        end
    end
    return res

end


function multi_dim_ryser(U, gram_matrix)

	@warn "obsolete, please convert to ryser_tensor(W::Array)"

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


function glynn(U; niter)

	"""approximates the permanent of U up to an additive error through the Gurvits/Glynn algorithm"""
	# see https://arxiv.org/abs/1212.0025

    function glynn_estimator(U,x)

        function product_U_x(U, x)
            result_product = one(typeof(U[1,1]))

            for j = 1:n
                result_product *= sum(U[j, :] .* x)
            end

            result_product
        end

        prod(x) * product_U_x(U, x)

    end

    n = size(U,1)
    result = zero(eltype(U))

    for i = 1:niter

        x = rand([-1,1], n)
        result += 1/niter * glynn_estimator(U,x)

    end

    result

end

function combine_glynn_estimators(est1, est2, n1, n2)

	"""combines two glynn estimations of n1,n2 trials to give one with (n1 + n2)"""

	(n1 * est1 + n2 * est2)/(n1 + n2)

end

function glynn_precision(U; rtol = 1e-5, miniter = 10^2, maxiter = 10^5, steps = 5)

	growth_factor = (maxiter/miniter)^(1/steps)

	growth_factor <= 2 ? (@warn "small growth factor, results may be inaccurate decrease steps") : nothing
	# the number of iterations is multiplied by this number at
	# each trial

	if miniter<1e2
		@warn "small miniter may give false results if the first estimations are by chance close to the next ones"
	end

	total_iter = 0
	estimates = []
	niters = []

	niter = miniter
	first_estimate = glynn(U; niter = niter)
	push!(estimates, first_estimate)
	push!(niters, niter)

	estimates[end]

	while total_iter < maxiter

		niter *= growth_factor
		new_estimate = glynn(U; niter = niter)
		push!(estimates, new_estimate)
		push!(niters, niter)
		rel_err = abs(estimates[end] - estimates[end-1]) / abs(0.5*(estimates[end] + estimates[end-1]))

		total_iter += niter

		if rel_err < rtol
			break
		end
	end

	estimates[end]

end

function number_of_steps(permanent_function, n)

    """gives the number of steps in an exact permanent function"""

    steps = Dict(ryser => n -> n*2^n, naive => n -> factorial(big(n)), naive_tensor => n -> factorial(big(n))^2, ryser_tensor => n -> 2^(2n))

    if permanent_function in keys(steps)
        steps[permanent_function](big(n))
    else
        error("not implemented")
    end
end

function error_permanent(permanent_function, n, T::Number)

    """gives an expected error on the permanent_function computed of a type T
    with a dimension n (matrix n*n, 3-tensor n*n*n)"""

    eps(T) * sqrt(number_of_steps(permanent_function, n))

end


end
