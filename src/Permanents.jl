module Permanents

export naive, naive_tensor, ryser, multi_dim_ryser, fast_glynn_perm, glynn, glynn_precision, ryser_tensor, positive_permanent

using Combinatorics
using LinearAlgebra
using ArgCheck
using Distributions

"""
Computes the permanent of the matrix ``U`` of dimension ``n`` using the definition
```math
perm(U) = ∑_{σ∈S_n}∏_{i=1}^n U_{i,σ(i)}
```
as a naive implementation in ``n!`` arithmetic operations.
"""
function naive(U::AbstractMatrix)

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

"""
Compute the permanent the permanent of a ``n^3``-dimensional 3-tensor of the form
```math
W_{k,l,j} = M_{k,j} M_{l,j}^* S_{l,k}
```
following Ryser's algorithm in approximatively ``2^{2n-1}`` iterations following
[Sampling of partially distinguishable bosons and the relation to the multidimensional permanent](https://arxiv.org/pdf/1410.7687.pdf).
!!! warning
	The current implementation does not use Gray code.
"""
function ryser_tensor(W::Array)

	@argcheck is_a_square_three_tensor(W) "tensor permanent implemented only for square 3-indices tensors"

    n = size(W)[1]
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

"""
Compute the permanent of a matrix ``A`` of dimension ``n`` using Ryser algorithm
with Gray ordering
```math
perm(A) = (-1)^n ∑_{S ⊆ 1 … n} (-1)^{|S|} ∏_{i=1}^n ∑_{j ∈ S} A_{i,j}
```
and time complexity ``O(2^{n-1}n)``.
!!! note
	source: [https://discourse.julialang.org/t/matrix-permanent/10766](https://discourse.julialang.org/t/matrix-permanent/10766)
"""
function ryser(A::AbstractMatrix)
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

"""
Compute the permanent of a matrix ``A`` of dimension ``n`` using Glynn formula
```math
perm(A) = (2^{n-1})^{-1} ∑_δ (∏_{k=1}^n δ_k) ∏_{j=1}^n ∑_{j=1}^n δ_i A_{i,j}
```
where ``δ ∈ (-1,+1)^n`` with time complexity ``O(n2^n)``.
!!! note
	source: [https://codegolf.stackexchange.com/questions/97060/calculate-the-permanent-as-quickly-as-possible](https://codegolf.stackexchange.com/questions/97060/calculate-the-permanent-as-quickly-as-possible)
"""
function fast_glynn_perm(U::AbstractMatrix{T}) where T

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

function postive_permanent(A::AbstractMatrix, niter=1e4)

	if !all(a -> a >= 0 for a in A)
		throw(ArgumentError("Input matrix must be with non-negative entries"))
	elseif size(A)[1] != size(A)[2]
		throw(ArgumentError("Input matrix must be a square matrix"))
	end

	n = size(A)[1]
	nb_success = 0

	function rescale_mat(A::AbstractMatrix, eps::Real)

		B = A
		n = size(A)[1]
		x = ones(n)
		y = x'

		row_sum = [sum(A[i,:]) for i in 1:n]

		while maximum([abs(r-1) for r in row_sum]) > eps
			x = x .* row_sum.^(-1)
			B = diagm(row_sum.^(-1)) * B
			col_sum = [sum(B[:,i]) for i in 1:n]
			y = y .* col_sum.^(-1)
			y = y[:,1]
			B *= diagm(col_sum.^(-1))
			row_sum = [sum(B[i,:]) for i in 1:n]
		end

		return [B, x, y]

	end

	function hl_factor(x::AbstractVector)

		x0 = [Bool(i==0) for i in x]
		x = x .+ 0.5*x0

		x_1 = [Bool(i>1) for i in x]
		x__1 = [Bool(i<=1) for i in x]
		log_x = [log(abs(i)) for i in x]

		return x_1 .* (x .+ 0.5.*log_x .+ exp(1) .- 1) + x__1 .* (1 .+ (exp(1)-1).*x)

	end

	B, x, y = rescale_mat(A, 1e-5)
	row_scaled = [maximum(B[i,:]) for i in 1:n].^(-1)
	C = diagm(row_scaled) * B
	C_i = copy(C)

	for i in 1:niter

		col = 1
		C = copy(C_i)
		row_sum = [sum(C[k,:]) for k in 1:n]

		while col <= n

			hl_i = hl_factor(row_sum)
			hl1 = prod(hl_i/exp(1))
			hl_i = hl_factor(row_sum - C[:,col])
		 	hl2 = prod(hl_i/exp(1))

			row_prob = exp(1) * C[:,col].*hl2/hl1./hl_factor(row_sum - C[:,col])
			val = rand()
			inter_bool = [Bool(sum_k < val) for sum_k in cumsum(row_prob)]
			row_samp = sum(inter_bool) + 1

			if row_samp == n+1
				col = n+2
				push!(row_sum, 0)
			else
				row_sum -= C[:,col]
				C[row_samp,:] = zeros(n)
				col += 1
				row_sum[row_samp] = 0
			end

		end

		if col == n+1
			nb_success += 1
		end

	end

	C = C_i
	row_sum = sum(C[:,k] for k in 1:n)
	hl_C = prod(hl_factor(row_sum)/exp(1))
	res = hl_C * nb_success/niter

	return res /prod(row_scaled)/prod(x)/prod(y)

end

function quantumInspired(A::AbstractMatrix, eps=1e-1, probFail=1e-1, C= 2)

	if !ishermitian(A)
		throw(ArgumentError("Input matrix must be hermitian"))
	elseif !isposdef(A)
		throw(ArgumentError("Input matrix must be at least positive semidefinite"))
	end

	D = eigvals(A)
	U = eigvecs(A)

	λ_max = maximum(D)
	C = 2
	n = length(D)

	rescaling_factor = (C * λ_max)^(2n) / prod((C*λ_max-λ) for λ in D)
	sample_size = trunc(Int, (rescaling_factor^2 * exp(-2n)) / (2eps^2) * log(1/probFail))
	weight_array = Vector{Number}(undef, sample_size)

	Threads.@threads for j in 1:sample_size

		sample_array = Vector{Number}(undef, n)

		for i in 1:n
			mean_i = D[i] / (1-D[i])
			d = Normal{Float64}(0.0, mean_i/2)
			sample_array[i] = rand(d) + rand(d)im
		end

		sampled_amplitudes = [dot(U[:,i], sample_array) for i in 1:n]
		weight_array[j] = prod(exp(-abs(sampled_amplitudes[i])^2 * abs(sampled_amplitudes[i])^2) for i = 1:n)

	end

	sample_mean = 1/sample_size * sum(weight_array[j] for j in 1:sample_size)
	return sample_mean * rescaling_factor

end


end
