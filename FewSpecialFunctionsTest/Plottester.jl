using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)


x = range(1,15,250)

plot(x,irregular_Coulomb.(1,0.5,x), label=L"G_1(0.5,ρ)")
plot!(x,irregular_Coulomb.(1,-0.5,x), label=L"G_1(-0.5,ρ)")
plot!(x,irregular_Coulomb.(1,0.3,x), label=L"G_1(0.3,ρ)")
plot!(x,irregular_Coulomb.(1,-0.3,x), label=L"G_1(-0.3,ρ)")
title!("Irregular Coulomb wave functions")
xlabel!(L"\rho")

