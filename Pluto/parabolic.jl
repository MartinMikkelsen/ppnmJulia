using FewSpecialFunctions, Plots, SpecialFunctions

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5,
    framestyle=:box,
    label=nothing,
    grid=true,
    palette=:tab10,
)

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

function confluent_hypergeometric_1F1_deriv(a,b,z)
    return a/b*confluent_hypergeometric_1F1(a+1,b+1,z)
end


x = range(0,10,100)
plot(x,confluent_hypergeometric_1F1_deriv.(1,3,x))

function f1(a,z)
    return exp(-0.25*z^2)*confluent_hypergeometric_1F1(0.5*a+0.25,0.5,0.5*z^2)
end

function f2(a,z)
    return z*exp(-0.25*z^2)*confluent_hypergeometric_1F1(0.5*a+0.75,1.5,0.5*z^2)
end

function parabolic_U(a,z)
    ξ = 0.5a+0.25
    return 1/(2^ξ*sqrt(π))*(cos(ξ*π)*gamma(0.5-ξ)*f1(a,z)-sqrt(2)*sin(ξ*π)*gamma(1-ξ)*f2(a,z))
end

function parabolic_V(a,z)
    ξ = 0.5a+0.25
    return 1/(2^ξ*sqrt(π)*gamma(0.5*-a))*(sin(ξ*π)*gamma(0.5-ξ)*f1(a,z)+sqrt(2)*cos(ξ*π)*gamma(1-ξ)*f2(a,z))
end

x = range(0.0,10,100)

plot(x,parabolic_U.(0.1,x))