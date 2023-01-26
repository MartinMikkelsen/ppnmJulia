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


x = range(0,25,1000)

plot(x,Debye_function.(1,x),label=L"D_1(x)")
plot!(x,Debye_function.(2,x),label=L"D_2(x)")
plot!(x,Debye_function.(3,x), label=L"D_3(x)")
title!("Debye Functions")
xlabel!(L"x")
savefig("./FewSpecialFunctionsTest/Debye_Example.pdf")