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


x = range(0,15,1000)

plot(x,Clausen.(x), label=L"Cl_2(ϕ)")
xlabel!(L"ϕ")
title!("Clausen function")
savefig("./FewSpecialFunctionsTest/Clausen_Example.pdf")

