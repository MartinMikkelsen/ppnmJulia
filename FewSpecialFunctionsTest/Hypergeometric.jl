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


function confluent_hypergeometric(a::Float64,b::Float64,z)
    if real(b)>real(a)>0
        M = 1/(gamma(a)*gamma(b-a))*quadgk(t -> exp(z*t)*t^(a-1)*(1-t)^(b-a-1),epsilon,1)[1]
        return M
    if b-a != 1,2,3
        M = gamma(1+a-b)/(2*π*1im*gamma(a))*quadgk(t -> exp(z*t)*t^(a-1)*(t-1)^(b-a-1),0,1)[1]
        return M
    if a != 1,2,3 > 0
        M = exp()
        return M
end

function hypergeometric(a,b,c,z)
    F = 1/(gamma(b)*gamma(c-b))*quadgk(t -> (t^(b-1)*(1-t)^(c-b-1))/((1-z*t)^(a)),1e-5,1)[1]
    return real(F)
end

x = range(0,10,100)

plot(x,confluent_hypergeometric.(1.6,3,x))