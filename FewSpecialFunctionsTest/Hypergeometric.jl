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

function hypergeometric_F(a,b,c,z)
    if (c-a-b) <= -1 
        return println("Divergence when (c-b-a) is smaller than or equal to -1. (c-b-a) = ", (c-b-a))
    end
    if (c-b-a) > 0 
        k = 1/(gamma(a)*gamma(c-b))
        if a == 1 && b == 1 && c == 2
            return -z^(-1)*log(1-z)
        end
        if b == 0.5+a && c == 2*a
            return 2^(2*a-1)*(1-z)^(-0.5)*(1+(1-z)^(0.5))^(1-2*a)
        end
        if z == 1 && (c-a-b) > 0 
            return (gamma(c)*gamma(c-a-b))/(gamma(c-a)*gamma(c-b))
        end
        if c == a-b+1 && z == -1
            return gamma(1+a-b)/(gamma(1+0.5*a-b)*gamma(0.5+0.5*a))
        end
        if a == 0.5 && b == 0.5 && c == 1.5 && z == z^2 
            return z^(-1)*asin(z)
        end
        if c>b>0 && abs(angle(1-z)) < π
            return k.*quadgk(t -> (t.^(b-1).*(1-t).^(c-b-1))/((1-z.*t).^(a)),1e-5,1.0)[1]
        end
    end
end

function Hypergeometric_0F1(b,z)
    if b == 1
        return 1/π * quadgk(t -> exp(2*sqrt(z)*cos(t)),0,π)[1]
    end
    if real(b)>0.5
        return (2*gamma(b))/(sqrt(π)*gamma(b-0.5)).*quadgk(t -> (1-t^2)^(b-1.5)*cosh(2*sqrt(z)*t),0,1)[1]
    end 
    if n typeof(val)<:Number > 0 
        return factorial(n-1)/(π)*z^((1-n)/2)*quadgk(t -> exp(2*sqrt(2)*cos(t))*cos((t-1)*t),0,π)[1]
    end

end

function confluent_hypergeometric_1F1(a,b,z)
    if real(a)>0
        return 1/(gamma(a)).*quadgk(t -> exp(-t)*t^(a-1).*Hypergeometric_0F1(b,z),0,Inf)[1]
    end
    if real(b)>real(a)>0
        return gamma(b)/(gamma(a).*gamma(b-a)).*quadgk(t -> exp(z*t).*t^(a-1).*(1-t)^(-a+b-1),0,1)[1]
    end
end

x = range(0,15,100)

plot(x,confluent_hypergeometric_1F1.(1.0,2,x))