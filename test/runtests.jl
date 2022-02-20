using Permanents
using Test

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

    U = fourier_matrix(7,normalized = false)
    @test glynn_precision(U; rtol = 1e-3, maxiter = 10^5) ≈ theoretical_permanents_fourier_matrix[7] rtol = 1e-2

end
