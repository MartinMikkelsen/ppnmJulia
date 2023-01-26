using Plots, FewSpecialFunctions, SpecialFunctions, QuadGK, BenchmarkTools, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)


x = range(-5,5,1000)

plot(x,Struve.(0,x),label=L"H_0(x)")
plot!(x,Struve.(1,x),label=L"H_1(x)")
plot!(x,Struve.(2,x),label=L"H_2(x)")
plot!(x,Struve.(3,x),label=L"H_3(x)")
plot!(x,Struve.(4,x),label=L"H_4(x)")
plot!(x,Struve.(5,x),label=L"H_5(x)")

xlabel!(L"x")
title!("Struve Functions")
