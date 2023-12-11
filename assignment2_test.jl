using Polynomials
using Random

# The actual test seed to be used is hidden
Random.seed!(99991099910101010)
rel_tol = 1e-6

include("assignment2_handout.jl")

# Test your implementation of Newton's form and Horner's method
function test_newton_int()
    N_test = 100
    n = rand(5:10)
    x = range(-1.0, stop=1.0, length=n)
    y = rand(length(x))

    # Use the fit function from Polynomials as ground truth
    p_true = fit(x, y)
    c_n = newton_int(x, y)
    X_eval = rand(N_test) * 2.0 .- 1.0
    vals = horner(c_n, x, X_eval)
    vals_true = p_true.(X_eval)

    # Count the number of points within tolerance
    score = sum(abs.((vals .- vals_true) ./ vals_true) .<= rel_tol)

    return score

end

score1 = test_newton_int()

# Test your implementation of Horner's method against a naive evaluation

function naive_polynomial_evaluation(c, x, X)
    n = length(c)
    m = length(X)
    p = c[1] * ones(m)
    for i in 1:m
        for j in 2:n
            p[i] += c[j] * prod([X[i] - x[k] for k in 1:j-1])
        end
    end
    return p
end

function test_horner()
    N_test = 100
    n = rand(5:10)
    x = range(-1.0, stop=1.0, length=n)
    y = rand(length(x))

    c_n = newton_int(x, y)
    X_eval = rand(N_test) * 2.0 .- 1.0
    vals_horner = horner(c_n, x, X_eval)
    vals_naive = naive_polynomial_evaluation(c_n, x, X_eval)
    # Count the number of points within tolerance
    score = sum(abs.((vals_horner .- vals_naive) ./ vals_naive) .<= rel_tol)

    return score
end

score2 = test_horner()

function time_horner()
    N_test = 100
    n = 55
    x = range(-1.0, stop=1.0, length=n)
    y = rand(length(x))

    c_n = newton_int(x, y)
    X_eval = rand(N_test) * 2.0 .- 1.0
    # Time both approaches
    t_horner = @elapsed horner(c_n, x, X_eval)
    println("t_horner: $(t_horner)")
    t_naive = @elapsed naive_polynomial_evaluation(c_n, x, X_eval)
    println("t_naive: $(t_naive)")
    score = 100 * (t_horner < 0.2 * t_naive)  # Yours should be at least 5 times as fast
    return score
end

score3 = time_horner()

println("Your total score on these test cases: $(score1+score2+score3)/300")
