using LinearAlgebra

""" 
Computes the coefficients of Newton's interpolating polynomial. 
    Inputs 
        x: vector with distinct elements x[i] 
        y: vector of the same size as x 
    Output 
        c: vector with the coefficients of the polynomial
"""
function newton_int(x, y)

    n = length(x)
    m = zeros(n, n)

    m[:, 1] .= y # Fill first column with y values

    # For each col, fill in values from i to n
    for i in 2:n # Col
        for j in i:n # Row
            m[j, i] = (m[j, i-1] - m[j-1, i-1]) / (x[j] - x[j-i+1])
        end
    end

    return diag(m)

end

"""
Evaluates a polynomial with Newton coefficients c 
defined over nodes x using Horner's rule on the points in X.
Inputs 
    c: vector with n coefficients 
    x: vector of n distinct points used to compute c in newton_int 
    X: vector of m points 
Output 
    p: vector of m points
"""
function horner(c, x, X)

    n, m = length(x), length(X)
    p = zeros(m)

    for i in 1:m

        p[i] = c[n]

        # O(n) time complexity
        for j in n-1:-1:1
            p[i] = c[j] + (X[i] - x[j]) * p[i]
        end
    end

    return p
end

"""
Computes the number of equally spaced points to use for 
interpolating cos(ω*x) on interval [a, b] for an absolute
error tolerance of tol.

Inputs
    a: lower boundary of the interpolation interval
    b: upper boundary of the interpolation interval
    ω: frequency of cos(ω*x)
    tol: maximum absolute error 
Output
    n: number of equally spaced points to use 	 
"""
function subdivide(a, b, ω, tol)
    for n in 2:100
        h = (b - a) / n

        curr = ω^(n + 1) * h^(n + 1) / (4 * (n + 1))

        if curr < tol
            # n is number of intervals
            # n + 1 is number of nodes
            return n + 1
        end
        n += 1
    end
end

"""
Computes Chebyshev nodes in the interval [a, b] for the function
cos(ω*x) for a maximum absolute error of tol.

Inputs
    a: lower boundary of the interpolation interval
    b: upper boundary of the interpolation interval
    ω: frequency of cos(ω*x)
    tol: maximum absolute error
Output
    x: distinct Chebyshev nodes	 
"""
function chebyshev_nodes(a, b, ω, tol)

    n = 2

    # Max value for factorial is 20
    while n < 20
        curr = ω^(n + 1) / 2^n / factorial(n + 1) * ((b - a) / 2)^(n + 1)

        if curr < tol
            break
        end
        n += 1
    end

    # n is the number of intervals
    # n + 1 is the number of nodes

    x = zeros(n + 1)

    for i in 1:n+1
        x[i] = 0.5 * (a + b) + 0.5 * (b - a) * cos(π * (2 * i - 1) / (2 * n + 2))
    end

    return x
end
