using Combinatorics
using LinearAlgebra
using ArgCheck
using Distributions
using StatsBase
using BenchmarkTools  

" To use 8 thread start JULIA as follows:
 julia --threads 8  
 https://docs.julialang.org/en/v1/manual/multi-threading/
"

function ryser_parallel_MultiThread(A::Matrix)

    "Computes the permanent of a general matrix n*n A 
	following the Ryser's algorithm by using Multi-Threading
    https://arxiv.org/pdf/1904.06229.pdf : Appendix A
    "
    n=size(A,1)
    m=min(n,Threads.nthreads())
    S=zeros(m)
    
    Threads.@threads for  k in 0:m-1
        S[k+1]=2*(-1)^n*sub_permanent(A,n,m,k)
        
    end
    
    return sum(S)
end




function ryser_parallel(A::Matrix)
    # Same function as "ryser_parallel" without the Multi-Threading
    n=size(A,1)
    m=min(size(A,1),Threads.nthreads())
    
    return sum([2*(-1)^n*sub_permanent(A,n,m,k) for  k in 0:m-1])
end



function sub_permanent(A::Matrix,n::Int,m::Int,k::Int)
    
    r,l=distribute(2^(n-1),m,k)
    x= (r ⊻ (r>>>1))
    x=decimal_to_binary_vector(x,2^(n-1))
    t=(-1)^(r)
    S=0
    w=[A[i,n]-(1/2)sum([A[i,j] for j in 1:n]) for i in 1:n]
    
    for j in 1:n-1
        if x[j]==1
            w=[w[i]+A[i,j] for i in 1:n]
        end
    end
    
    for e in 1:l
        p=prod(w)
        t,j,x=next_setp(t,x)
        S+=t*p
        z=2*x[j]-1
        w=[w[i]+z*A[i,j] for i in 1:n]
    end
    
    return S
end

function distribute(n::Int,m::Int,i::Int)

    q=n÷m
    r=n-q*m
    k=(i)*q+min(i,r)
    (i<r) ? l=q+1 : l=q
    return k,l 
end

function next_setp(t,x)
    j=1
    t=-t
    
    if t==1
        while x[j]==0
            j=j+1
        end
        j=j+1
    end
    x[j]=1-x[j]
   
    return t,j,x
end

function decimal_to_binary_vector(n::Int,m::Int)

    i = 1;
    bn=[0 for j in 1:floor(Int,log2(m))+1]
    while (n > 0)
        bn[i] = n % 2
        n=n÷2
        i=i+1
        
    end
    bn

end

