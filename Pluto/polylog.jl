using Plots
using SpecialFunctions

function binomial_coefficient(n, j)
    return binomial(BigInt(n), BigInt(j))
end

function h(s, z, n)
    return (z^2 / (4 * (z - 1)))^n * z / gamma(s)
end

function find_n_for_error_bound(z, s, epsilon)
    n = 1
    while true
        error_bound = h(s, z, n)
        if error_bound < epsilon
            break
        end
        n += 1
    end
    return n
end

function f(z, s, n)
    result = 0.0
    
    # First sum
    for k = 1:n
        result += z^k / k^s
    end
    
    # Second sum with nested loops
    second_sum = 0.0
    for k = n+1:2n
        inner_sum = 0.0
        for j = 0:(2n - k)
            inner_sum += (-z)^j * binomial_coefficient(2n - k, j)
        end
        second_sum += z^k / k^s * inner_sum
    end
    
    result += 1 / (1 - z)^n * second_sum + h(s, z, n)
    
    return result
end

# Define parameters
s = 2.0
epsilon = 1e-3
z_values = -4.999:0.001:0.99  # Exclude the problematic region

# Calculate f(z, s) for each z value
f_values = [f(z, s, find_n_for_error_bound(z, s, epsilon)) for z in z_values]

# Create and show the plot
plot(z_values, f_values, xlabel="z", ylabel="f(z, s)", label="f(z, s)", legend=true)
