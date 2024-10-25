using Permanents
using Test
using LinearAlgebra
using DynamicPolynomials

const ATOL = 1000eps()

function fourier_matrix(n; normalized = true)
    U = Array{ComplexF64}(undef, n, n)
    for i in 1:n
        for j in 1:n
            U[i,j] = exp((2im * pi / n) * (i-1) * (j-1))
        end
    end
    if normalized
        return 1/sqrt(n) * U
    else
        return U
    end
end

@testset "Permanents.jl" begin
    # Known theoretical values for Fourier matrix permanents
    theoretical_permanents_fourier_matrix = [1,0,-3,0,-5,0,-105,0,81,0,6765,0,175747,0,30375,0,25219857,0,142901109,0,4548104883,0]

    @testset "Basic permanent calculations" begin
        for n = 1:6
            @test naive(ones(Int,n,n)) == factorial(n)
            @test naive_tensor(ones(n,n,n)) ≈ factorial(n)^2 atol=ATOL
            
            U = fourier_matrix(n, normalized = false)
            @test naive(U) ≈ theoretical_permanents_fourier_matrix[n] atol=ATOL
            @test ryser(U) ≈ theoretical_permanents_fourier_matrix[n] atol=ATOL
            @test fast_glynn_perm(U) ≈ theoretical_permanents_fourier_matrix[n] atol=ATOL
            @test multi_dim_ryser(U, ones(n,n)) ≈ abs(theoretical_permanents_fourier_matrix[n])^2 atol=ATOL
        end
    end

    @testset "Incomplete Rank Algorithm" begin
        @testset "Standard cases" begin
            test_sizes = [(4,3), (6,4), (5,4), (5,2)]
            
            for (n, r) in test_sizes
                @testset "$(n)×$(n) matrix, rank $r" begin
                    for _ in 1:5
                        A = rand_gram_matrix_rank(n, r)
                        perm_ryser = ryser(A)
                        perm_incomplete = incomplete_rank(A)
                        @test perm_ryser ≈ perm_incomplete atol=ATOL
                    end
                end
            end
        end

        @testset "Special matrices" begin
            # Zero matrix
            A = zeros(ComplexF64, 4, 4)
            @test ryser(A) ≈ incomplete_rank(A) atol=ATOL

            # Rank deficient diagonal matrix
            A = Diagonal([1.0, 1.0, 1.0, 0.0])
            @test ryser(A) ≈ incomplete_rank(A) atol=ATOL

            # Complex matrix
            A = rand_gram_matrix_rank(4, 2) + im * rand_gram_matrix_rank(4, 2)
            @test ryser(A) ≈ incomplete_rank(A) atol=ATOL
        end

        @testset "Error cases" begin
            # Non-square matrix
            A = randn(3, 4)
            @test_throws ArgumentError ryser(A)
            @test_throws DimensionMismatch incomplete_rank(A)

            # Full rank warning
            A = Matrix(1.0I, 4, 4)
            @test_logs (:warn, "while the algorithm will work for any rank, you should rather use another one if full rank, using ryser instead") incomplete_rank(A)
        end

        @testset "Fourier matrix comparisons" begin
            for n = 2:4
                U = fourier_matrix(n, normalized=false)
                # Force rank deficiency by zeroing last column
                U[:, end] .= 0
                @test ryser(U) ≈ incomplete_rank(U) atol=ATOL
            end
        end
    end
end