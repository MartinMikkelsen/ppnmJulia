using Plots, FewSpecialFunctions, SpecialFunctions, QuadGK, BenchmarkTools, LaTeXStrings, PyPlot

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

function H0(x::Float64)
    if abs(x) < 0.1
        # series expansion for small x
        h0 = x^2/3 - x^4/30 + x^6/840
        # add more terms to the series expansion if needed
        return h0
    else
        return besselj(0, x) * cos(x) - bessely(0, x) * sin(x)
    end
end

function H1(x::Float64)
    if abs(x) < 0.1
        # series expansion for small x
        h1 = x - x^3/6 + x^5/60 - x^7/2520
        # add more terms to the series expansion if needed
        return h1
    else
        return besselj(1, x) * cos(x) - bessely(1, x) * sin(x)
    end
end

x = range(-5,5,1000)

plot(x,Struve.(0,x),label=L"H_0(x)")
plot!(x,Struve.(1,x),label=L"H_1(x)")
plot!(x,Struve.(2,x),label=L"H_2(x)")
plot!(x,Struve.(3,x),label=L"H_3(x)")
plot!(x,Struve.(4,x),label=L"H_4(x)")
plot!(x,Struve.(5,x),label=L"H_5(x)")

xlabel!(L"x")
title!("Struve Functions")
