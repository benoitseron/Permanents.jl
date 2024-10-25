using ArgCheck
using LinearAlgebra
using DynamicPolynomials
using Combinatorics  # Required for factorial calculations
using MultivariatePolynomials
using Permanents

# Function to generate a random Gram matrix with specified rank
function rand_gram_matrix_rank(n, r)
    function normalized_random_vector(r)
        v = rand(ComplexF64, r)
        1/norm(v) .* v
    end

    generating_vectors = hcat([normalized_random_vector(r) for i = 1:n]...)
    generating_vectors' * generating_vectors
end



# Example usage

A = rand_gram_matrix_rank(4, 3)
permanent_A = ryser(A)
println("Permanent of matrix A: ", permanent_A)

incomplete_rank(A)

 # Test case 1: Small matrix (4×4) with rank 3
 @testset "4×4 matrix, rank 3" begin
    for _ in 1:5  # Run multiple times with different random matrices
        A = rand_gram_matrix_rank(4, 3)
        perm_ryser = ryser(A)
        perm_incomplete = incomplete_rank(A)
        @test isapprox(perm_ryser, perm_incomplete, rtol=1e-10)
    end
end




@testset "Permanent Calculations - Incomplete Rank vs Ryser" begin
        
    # Test case 1: Small matrix (4×4) with rank 3
    @testset "4×4 matrix, rank 3" begin
        for _ in 1:5  # Run multiple times with different random matrices
            A = rand_gram_matrix_rank(4, 3)
            perm_ryser = ryser(A)
            perm_incomplete = incomplete_rank(A)
            @test isapprox(perm_ryser, perm_incomplete, rtol=1e-10)
        end
    end

    # Test case 2: Larger matrix (6×6) with rank 4
    @testset "6×6 matrix, rank 4" begin
        for _ in 1:5
            A = rand_gram_matrix_rank(6, 4)
            perm_ryser = ryser(A)
            perm_incomplete = incomplete_rank(A)
            @test isapprox(perm_ryser, perm_incomplete, rtol=1e-10)
        end
    end

    # Test case 3: Edge case - Almost full rank
    @testset "5×5 matrix, rank 4" begin
        for _ in 1:5
            A = rand_gram_matrix_rank(5, 4)
            perm_ryser = ryser(A)
            perm_incomplete = incomplete_rank(A)
            @test isapprox(perm_ryser, perm_incomplete, rtol=1e-10)
        end
    end

    # Test case 4: Low rank case
    @testset "5×5 matrix, rank 2" begin
        for _ in 1:5
            A = rand_gram_matrix_rank(5, 2)
            perm_ryser = ryser(A)
            perm_incomplete = incomplete_rank(A)
            @test isapprox(perm_ryser, perm_incomplete, rtol=1e-10)
        end
    end

    # Test case 5: Special matrices
    @testset "Special matrices" begin
        # Zero matrix
        A = zeros(ComplexF64, 4, 4)
        @test isapprox(ryser(A), incomplete_rank(A), rtol=1e-10)
        
        # Diagonal matrix with rank deficiency
        A = Diagonal([1.0, 1.0, 1.0, 0.0])
        @test isapprox(ryser(A), incomplete_rank(A), rtol=1e-10)
        
        # Complex matrix
        A = rand_gram_matrix_rank(4, 2) + im * rand_gram_matrix_rank(4, 2)
        @test isapprox(ryser(A), incomplete_rank(A), rtol=1e-10)
    end

    # Test case 6: Error handling
    @testset "Error handling" begin
        # Non-square matrix
        A = randn(3, 4)
        @test_throws ArgumentError ryser(A)
        @test_throws DimensionMismatch incomplete_rank(A)
        
        # Full rank matrix should issue warning
        A = Matrix(1.0I, 4, 4)
        @test_logs (:warn, "while the algorithm will work for any rank, you should rather use another one if full rank, using ryser instead") incomplete_rank(A)
    end
end