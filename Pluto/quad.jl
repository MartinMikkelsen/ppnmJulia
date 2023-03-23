function transform_limit(f, a)
    return t -> f(a + t/(t-1)*1/(1-t)^2)
end

function recursive_integrator(f, a, b=Inf, atol=1e-8, rtol=1e-6)
    # Use an adaptive trapezoidal rule to approximate the integral
    function trap_rule(f, a, b)
        return (b-a)*(f(a) + f(b))/2
    end
    
    # Apply a transformation to handle infinite limits
    if a == -Inf && b == Inf
        f, a, b = transform_limit(f, 0)
    elseif a == -Inf
        f, a, b = transform_limit(f, b)
    elseif b == Inf
        f, a, b = transform_limit(f, a)
    end

    # Recursive function to subdivide the interval and estimate the error
    function recurse(f, a, b, fa, fb, fc, tol)
        c = (a+b)/2
        fd = f(c)
        fab = trap_rule(f, a, b)
        fac = trap_rule(f, a, c)
        fcb = trap_rule(f, c, b)
        err = abs(fab - (fac + fcb))
        if err <= tol
            return fac + fcb
        else
            result1 = recurse(f, a, c, fa, fd, fc, tol/sqrt(2))
            result2 = recurse(f, c, b, fc, fd, fb, tol/sqrt(2))
            return result1 + result2
        end
    end
    
    # Evaluate the integral and return the result
    fa = f(a)
    fb = f(b)
    fc = f((a+b)/2)
    result = recurse(f, a, b, fa, fb, fc, atol + rtol*abs(fb))
    return result
end

# Define a function with infinite limits of integration
f(x) = 1/x^2

# Compute the integral from 0 to infinity using the exponential transform
result = recursive_integrator(f, 1.1)

# Print the result
println("Result: ", result)
