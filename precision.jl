function number_of_steps(permanent_function, n)

    """gives the number of steps in an exact permanent function"""

    if permanent_function == ryser
        n*2^n
    elseif permanent_function == naive
        factorial(big(n))
    elseif permanent_function == naive_tensor
        factorial(big(n))^2
    else
        error("not implemented")
    end
end
