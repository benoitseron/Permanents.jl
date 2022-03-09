using Permanents
using Test
using LinearAlgebra

const ATOL = 1e6*eps()

function fourier_matrix(n; normalized = true)

    U = Array{ComplexF64}(undef, n, n)

    for i in 1:n
        for j in 1:n
            U[i,j] = exp((2im * pi / n)* (i-1) * (j-1))
        end
    end

    if normalized
        return 1/sqrt(n) * U
    else
        return U
    end
end

function rand_haar(n::Int)

    """generates a Haar distributed unitary matrix of size n*n"""
    # follows https://case.edu/artsci/math/esmeckes/Meckes_SAMSI_Lecture2.pdf

    qr(randn(ComplexF64, n,n)).Q
end


function rand_gram_matrix(n::Int)

	"""random gram matrix with full rank """

	function normalized_random_vector(n)
		v = rand(ComplexF64, n)
		1/norm(v) .* v
	end

	generating_vectors = 0.
	while det(generating_vectors) ≈ 0. #check it's a basis
	    generating_vectors = hcat([normalized_random_vector(n) for i = 1:n]...)

	end

	generating_vectors' * generating_vectors

end


@testset "Permanents.jl" begin
    theoretical_permanents_fourier_matrix = [1,0,-3,0,-5,0,-105,0,81,0,6765,0,175747,0,30375,0,25219857,0,142901109,0,4548104883,0]
    #from "on the permanent of schur matrix", graham

    for n = 1:6

        @test naive(ones(Int,n,n)) == factorial(n)
        @test naive_tensor(ones(n,n,n)) ≈ factorial(n)^2

        U = fourier_matrix(n, normalized = false)
        @test naive(U) ≈ theoretical_permanents_fourier_matrix[n] atol=ATOL
        @test ryser(U) ≈ theoretical_permanents_fourier_matrix[n] atol=ATOL
        @test fast_glynn_perm(U) ≈ theoretical_permanents_fourier_matrix[n] atol=ATOL
        @test multi_dim_ryser(U, ones(n,n)) ≈ abs(theoretical_permanents_fourier_matrix[n]).^2 atol=ATOL
    end

	### consistency of ryser_tensor(W) ≈  multi_dim_ryser(U, S) ###

    U = fourier_matrix(7,normalized = false)
    @test glynn_precision(U; rtol = 1e-3, maxiter = 10^5) ≈ theoretical_permanents_fourier_matrix[7] rtol = 1e-2

	n = 7
	U = rand_haar(n)
	S = rand_gram_matrix(n)

	W = Array{eltype(S)}(undef, (n,n,n))

	for ss in 1:n
	    for rr in 1:n
	        for j in 1:n
	            W[ss,rr,j] = U[ss, j] * conj(U[rr, j]) * S[rr, ss]
	        end
	    end
	end

	@test ryser_tensor(W) ≈  multi_dim_ryser(U, S)

	U = rand(4,4)
	@test positive_entry(U) ≈ naive(U) atol=1e-2

end
