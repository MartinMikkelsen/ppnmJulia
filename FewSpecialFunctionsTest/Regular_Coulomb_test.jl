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

plot(x,regular_coulomb.(0,0.3,x), label=L"F_0(0.3,ρ)")
plot!(x,regular_coulomb.(0,-0.3,x), label=L"F_0(0.3,ρ)")
xlabel!(L"ρ")
title!("Regular Coulomb Wave Functions")
savefig("./FewSpecialFunctionsTest/Regular_Coulomb1.pdf")

plot(x,regular_coulomb.(1e-5,5,x), label=L"F_0(5,ρ)")
plot!(x,regular_coulomb.(1,5,x), label=L"F_1(5,ρ)")
plot!(x,regular_coulomb.(2,5,x), label=L"F_2(5,ρ)")
plot!(x,regular_coulomb.(3,5,x), label=L"F_3(5,ρ)")
title!("Regular Coulomb Wave Functions")
xlabel!(L"ρ")
savefig("./FewSpecialFunctionsTest/Regular_Coulomb2.pdf")