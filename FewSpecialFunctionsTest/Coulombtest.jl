using FewSpecialFunctions, SpecialFunctions, Plots, QuadGK

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

function regular_coulomb(ℓ,η,ρ)
    First = ρ.^(ℓ+1)*2^ℓ.*exp(1im.*ρ-(π.*η/2))/(abs(gamma(ℓ+1+1im*η)))
    Integral_value = First.*complex_quadrature(t -> exp(-2*1im.*ρ*t).*t.^(ℓ+1im*η)*(1-t).^(ℓ-1im*η),1e-6,1)
    return Integral_value
end

#Needs work
function C(ℓ,η)
    return 2^ℓ*exp(-π*η/2).*(abs(gamma(ℓ+1+1im*η))/(factorial(2*ℓ+1)))
end

function irregular_coulomb(ℓ,η,ρ)
    First = exp(-1im*ρ)*ρ^(-ℓ)/(factorial(2*ℓ+1).*C(ℓ,η))
    if ℓ == 0
        Integral_value = complex_quadrature(t -> exp(-t)*t^(-ℓ-1im*η)*(t+2*1im*ρ)^(ℓ+1im*η),0.0,Inf)
        return First.*Integral_value
    end
end

function θ(ℓ,η,ρ)
    return ρ - η.*log(2*ρ)-0.5*ℓ*π+angle.(gamma(ℓ+1+1im*η))
end

function Coulomb_H(ℓ,η,ρ)
    return (exp(-1im*ρ)*ρ^(-ℓ))/(factorial(2*ℓ+1)*C(ℓ,η)).*complex_quadrature(t -> exp(-t)*t^(ℓ-1im*η)*(t+2*1im*ρ)^(ℓ+1im*η),0,Inf)
end 

x = range(1,15,500)
x2 = range(1,2.5,500)

ℓ, η = 0, 2

plot(x,real(Coulomb_H.(ℓ,η,x)), label="G")

plot!(x,exp.(-θ.(ℓ,η,x)),label="large ρ")
#plot!(x2,(x2.^(-ℓ))/((2*ℓ+1).*C(ℓ,η)), label="small ρ")

function regular_coulomb_1F1(ℓ,η,ρ)
    return C_test(ℓ,η).*ρ^(ℓ+1).*exp(1im*ρ).*confluent_hypergeometric_1F1(ℓ+1+1im*η,2*ℓ+2,-2*1im*ρ)
end