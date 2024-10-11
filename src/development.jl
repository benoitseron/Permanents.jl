using ArgCheck
using LinearAlgebra
using Polynomials

function rand_gram_matrix_rank(n, r)
    """Generates a random Gram matrix of size n×n with rank at most r."""
    function normalized_random_vector(r)
        v = rand(ComplexF64, r)
        v / norm(v)  # Normalize the vector
    end

    generating_vectors = hcat([normalized_random_vector(r) for _ in 1:n]...)
    return generating_vectors' * generating_vectors
end

function incomplete_rank(A::Matrix; atol = 1000eps())
    """Computes the permanent of a general n×n matrix A of rank r < n
    using the polynomial approach described in the algorithm."""
    
    n = size(A, 1)
    r = rank(A, atol = atol)

    @argcheck r < n  # Ensure that rank is less than n

    qr_dec = qr(A)
    G = qr_dec.Q
    B = qr_dec.R

    # Step 2: Define the polynomials L and R
    L_poly = prod([Polynomial(B[:, j][1:r]) for j in 1:n]) ### not good, it's not the same x!
    R_poly = prod([Polynomial(G[i, :][1:r]) for i in 1:n])

    # Expand the polynomials L and R into their coefficients
    L_coeffs = L_poly.coeffs
    R_coeffs = R_poly.coeffs

    # Step 3: Compute the permanent using the coefficients from the expanded polynomials
    permanent_value = sum(L_coeffs[k] * R_coeffs[k] for k in 1:length(L_coeffs)) 

	############# the above is missing the factorials!

    return permanent_value
end

# Example Usage
atol = 1000eps()
A = rand_gram_matrix_rank(3, 2)
incomplete_rank_value = incomplete_rank(A, atol=atol)
ryser_value = ryser(A)

incomplete_rank_value ≈ ryser_value || error("The computed permanents do not match.")

println("The computed permanent of the matrix is: ", permanent_value)

qr_dec = qr(A)
    G = qr_dec.Q
    B = qr_dec.R

	G * B

	A

B[1,3]


# Step 2: Define the polynomials L and R
n = 3
r = 2

L_poly = prod([Polynomial(B[:, j][1:r]) for j in 1:n])
R_poly = prod([Polynomial(G[i, :][1:r]) for i in 1:n])

L_coeffs = L_poly.coeffs
R_coeffs = R_poly.coeffs

j = 2

Polynomial(B[:, j][1:r])