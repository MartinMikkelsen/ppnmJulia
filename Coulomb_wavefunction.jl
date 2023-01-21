using QuadGK, SpecialFunctions, Plots, BenchmarkTools

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

function regular_coulomb(ℓ,η,ρ)
    First = ρ^(ℓ+1)*2^ℓ*exp(1im*ρ-(π*η/2))/(abs(gamma(ℓ+1+1im*η)))
    Integral_value = complex_quadrature(t -> exp(-2*1im.*ρ*t)*t^(ℓ+1im*η)*(1-t)^(ℓ-1im*η),0,1)
    return First.*Integral_value
end

function C(ℓ,η)
    return 2^ℓ*exp(-π*η/2)*(abs(gamma(ℓ+1+1im*η))/(factorial(2*ℓ+1)))
end

function irregular_coulomb(ℓ,η,ρ)
    First = exp(-1im*ρ)*ρ^(-ℓ)/(factorial(2*ℓ+1)*C(ℓ,η))
    Integral_value = complex_quadrature(t -> exp(-t)*t^(-ℓ-1im*η)*(t+2*1im*ρ)^(ℓ+1im*η),0.01,1e4)
    return First.*Integral_value
end


xes = range(0,10,100)
coulombwave1 = real([regular_coulomb(1,-2,i) for i in xes])
coulombwave2 = real([regular_coulomb(2,-10,i) for i in xes])


plot!(xes,coulombwave1, linewidth=2.5, label="Regular Coulomb wave(1,-2,ρ)")
plot!(xes,coulombwave2, linewidth=2.5, label="Regular Coulomb wave(2,-10,ρ)")


using Test

@testset "complex_quadrature" begin
    @test_approx_eq complex_quadrature(x->x^2,0,1) ≈ 0.3333333333333333im
    @test_approx_eq complex_quadrature(x->x^3,0,1) ≈ 0.25im
end
