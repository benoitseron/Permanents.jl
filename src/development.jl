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

