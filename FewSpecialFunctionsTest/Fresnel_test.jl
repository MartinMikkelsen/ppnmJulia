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


x = range(-25,25,1000)

plot(x,Fresnel_C_integral.(x),label=L"C(x)")
plot!(x,Fresnel_C_err.(x), ls=:dash, lw=1.5, label=L"\tilde{C}(x)")
title!("Fresnel Integral")
xlabel!(L"x")
#savefig("./FewSpecialFunctionsTest/Fresnel_C_example.pdf")

plot(x,Fresnel_S_integral.(x),label=L"S(x)")
plot!(x,Fresnel_S_err.(x), ls=:dash, lw=1.5, label=L"\tilde{S}(x)")
title!("Fresnel Integral")
xlabel!(L"x")
#savefig("./FewSpecialFunctionsTest/Fresnel_S_example.pdf")