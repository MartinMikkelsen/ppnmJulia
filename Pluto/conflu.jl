using FewSpecialFunctions, Plots, SpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5,
    framestyle=:box,
    label=nothing,
    grid=true,
    palette=:tab10,
)

@doc raw"""
    confluent_hypergeometric_1F1(a,b,z)
    
Returns the Kummer confluent hypergeometric function 

```math
    {}_1 F_1(a,b,z) = \sum_{k=0}^{\infty} \frac{(a)_k z^k}{(b)_k  k!}
```
"""
function confluent_hypergeometric_1F1(a,b,z)
    result = 1.0
    term = 1.0
    for n = 1:100
        term *= (a+n-1) / ((b+n-1)*n) * z
        result += term
        if abs(term) < 1e-12 # check if the term is negligible
            break
        end
    end
    return result
end

function confluent_hypergeometric_U(a, b, z)
    if b == a+1
        return z^(-a)
    else 
    f1 = confluent_hypergeometric_1F1(a, b, z)
    f2 = confluent_hypergeometric_1F1(1+a-b, 2-b, z)
    return gamma(1-b)/(gamma(1+a-b))*f1 +gamma(b-1)/(gamma(a))*f2
    end
end

println(confluent_hypergeometric_U(2.0,3.0,4.0))
#println(confluent_hypergeometric_U(0.25,5,1000))

