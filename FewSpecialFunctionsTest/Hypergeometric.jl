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

@doc raw"""
    hypergeometric_0F1(b,z)

Returns the confluent hypergeometric function given by

```math
    {}_0 F_1(a,b)
```
for the parameters ``a`` and ``b``
"""
function hypergeometric_0F1(b,z)
    if b == 1
        return 1/π .* quadgk(t -> exp(2*sqrt(z)*cos(t)),0,π)[1]
    end
    if real(b)>=0.5
        return (2*gamma(b))/(sqrt(π)*gamma(b-0.5)).*quadgk(t -> (1-t^2)^(b-1.5)*cosh(2*sqrt(z)*t),0,1)[1]
    end 
end
@doc raw"""
    confluent_hypergeometric_1F1(a,b,z)
    
Returns the Kummer confluent hypergeometric function 

```math
    {}_1 F_1
```
"""
function confluent_hypergeometric_1F1(a,b,z)
    if real(b)>real(a)>0
        return gamma(b)/(gamma(a)*gamma(b-a))*quadgk(t -> exp(z*t)*t^(a-1)*(1-t)^(-a+b-1),0,1)[1]
    end
    if real(a)>0 
        return 1/(gamma(a))*quadgk(u -> exp(-u/(1-u))*(u/(1-u))^(a-1)*hypergeometric_0F1(b,z*u/(1-u))*u/(1-u)^2,0,1)[1]
    end
end
@doc raw"""
    confluent_hypergeometric_U(a,b,z)
    
Returns the Kummer confluent hypergeometric function 

```math
    U(a,b,z)
```
"""
function confluent_hypergeometric_U(a,b,z)
    if real(z)>0 && real(a) > 0 
        return 1/gamma(a)*quadgk(u -> exp(-z*u/(1-u))*(u/(1-u))^(a-1)*(u/(1-u)+1)^(-a+b-1)*u/(1-u).^2,0,1)[1]
    end
end

function parabolic_cylinder_U(a,z)
    if real(a) > -0.5
        return exp(-0.25*z^2)/(gamma(0.5+a)).*quadgk(u -> (u/(1-u))^(a-0.5)*exp(-0.5*(u/(1-u))^2-z*u/(1-u))*1/(1-u)^2,0,1)[1]
    end
end

function parabolic_cylinder_D(ν,z)
    if real(ν) < 0 && real(z^2) > 0
        return 2^((ν+7)/2)/(gamma(-ν/2))*quadgk(u -> (u/(1-u))^(-ν/2-1)*(2*u/(1-u)+1)^((ν-1)/2)*exp(-z^2*u/(1-u)),0,1)[1]
    end
    if real(ν) > -1
        return sqrt(2)/(sqrt(π))*exp(z^2*0.25)*quadgk(u -> (u/(1-u))^ν*exp(-(u/(1-u))^2*0.5)*cos(z*u/(1-u)-π*ν*0.5),0,1)[1]
    end
end

function testing(a,b,z)
    return gamma(b-1)/(gamma(a)).*z.^(1-b).*confluent_hypergeometric_1F1(a-b+1,2-b,z)+gamma(1-b)/gamma(a-b+1).*confluent_hypergeometric_1F1(a,b,z)
end

x = range(0,15,500)

plot(x,hypergeometric_0F1.(1,x))
plot!(x,hypergeometric_0F1.(2,x))
plot!(x,hypergeometric_0F1.(3,x))