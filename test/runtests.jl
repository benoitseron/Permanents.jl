using Main.Permanents
using Test

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
    end

end
