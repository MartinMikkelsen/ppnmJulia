using Plots, FewSpecialFunctions, SpecialFunctions, QuadGK

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

function complex_quadrature(func,a,b)
    function real_func(x)
        return real(func(x))
    end
    function imag_func(x)
        imag(func(x))
    end
    real_integral = quadgk(real_func,a,b)
    imag_integral = quadgk(imag_func,a,b)
    return real_integral[1] + 1im*imag_integral[1]
end 

function confluent_hypergeometric_M(a,b,z)
    if a == 0
        return 1.0
    end 
    if a == 1 && b == 2
        return (exp(z)-1)/z
    end
    if a == 1 && b == 3
        return factorial(2)*(exp(z)-1-z)/(z^2)
    end
    if b==a
        return exp(z)
    end
    if b>a>0
        return 1/(gamma(a)*gamma(b-a)).*quadgk(t -> exp(z*t)*t^(a-1)*(1-t)^(b-a-1),0,1,rtol=1e-9)[1]
    end
    if b==a>0
        return exp(z)*quadgk(t^(-a)*exp(-t),z,Inf)[1]
    end
end

function confluent_hypergeometric_U(a,b,z)
    return 1/(gamma(a)).*quadgk(t -> exp(-z*t)*t^(a-1)*(1+t)^(b-a-1),0,Inf)[1]
end

function KummerU2(a,b,x)
tol = 1e-15
term = x*a/b
f = 1 + term
n = 1
an = a
bn = b
nmin = 10
    while(n < nmin) || max(abs(term) > tol)
        n = n + 1
        an = an + 1
        bn = bn + 1
        term = x.*term*an/bn/n
        f = f + term
    end
    return f
end