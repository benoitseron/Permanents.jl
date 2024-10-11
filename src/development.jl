using ArgCheck
using LinearAlgebra
using DynamicPolynomials
using Combinatorics  # Required for factorial calculations

# Function to generate a random Gram matrix with specified rank
function rand_gram_matrix_rank(n, r)
    function normalized_random_vector(r)
        v = rand(ComplexF64, r)
        1/norm(v) .* v
    end

    generating_vectors = hcat([normalized_random_vector(r) for i = 1:n]...)
    generating_vectors' * generating_vectors
end

# Function to compute the multinomial coefficient for a given tuple of integers
function multinomial_coeff(alpha::Vector{Int})
    return factorial(sum(alpha)) / prod(factorial(a) for a in alpha)
end

# Function to compute the permanent using Barvinok's approach
function compute_permanent(A::Matrix; atol = 1000eps())
    n = size(A, 1)
    r = rank(A, atol = atol)
    
    @argcheck r < n

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
    L_coeffs, L_monomials = coefficients(L_poly, x)
    R_coeffs, R_monomials = coefficients(R_poly, x)

    # Compute the permanent as sum of products of corresponding coefficients with multinomial weights
    perm_A = 0
    for i in 1:length(L_coeffs)
        alpha = exponents(L_monomials[i])  # Exponents of the monomial
        weight = multinomial_coeff(alpha)  # Multinomial coefficient
        perm_A += weight * L_coeffs[i] * R_coeffs[i]
    end

    return perm_A
end

# Example usage
A = rand_gram_matrix_rank(4, 3)
permanent_A = compute_permanent(A)
println("Permanent of matrix A: ", permanent_A)

atol = 1e-6

n = size(A, 1)
r = rank(A, atol = atol)

@argcheck r < n

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
L_coeffs, L_monomials = coefficients(L_poly, x)
R_coeffs, R_monomials = coefficients(R_poly, x)

# Compute the permanent as sum of products of corresponding coefficients with multinomial weights
perm_A = 0
for i in 1:length(L_coeffs)
    alpha = exponents(L_monomials[i])  # Exponents of the monomial
    weight = multinomial_coeff(alpha)  # Multinomial coefficient
    perm_A += weight * L_coeffs[i] * R_coeffs[i]
end

return perm_A

coefficients(L_poly)
L_poly

# L_poly = (0.0016178527654873668 - 0.0005447625684290993im)x₁x₂x₃² + (-0.007876692478264337 + 0.014851162367459429im)x₁x₂²x₃ + (-0.015110921413782563 - 0.04220793243383275im)x₁x₂³ + (0.0077351553445978135 - 0.003661296181176763im)x₁²x₃² + (0.07182552832945983 + 0.033234607197451374im)x₁²x₂x₃ + (-0.3173927489671689 + 0.30589499024755046im)x₁²x₂² + (0.48286900410494205 - 0.28568907287070455im)x₁³x₃ + (0.8153405008646424 + 1.7584828430813846im)x₁³x₂ + (7.442190392157916 - 5.390038134840203im)x₁⁴
#  need to extract each term and its coefficient