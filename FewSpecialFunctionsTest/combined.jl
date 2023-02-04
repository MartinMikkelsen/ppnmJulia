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

x = range(0.0,12.5,5000)
x2 = range(0.0,25,5000)
x3 = range(-25.0,25.0,5000)
x4 = range(0.0,12.5,5000)

y1 = Clausen.(x)
y2 = regular_coulomb.(1,0.3,x2)
y3 = Fresnel_S_err.(x3)
y4 = Struve.(1,x)

p1 = Plots.plot(x,y1, label=L"Cl_2(x)", title="Clausen")
p2 = Plots.plot(x2,y2, label=L"F_1(0.3,x)", title="Regular Coulomb")
p3 = Plots.plot(Fresnel_C_err.(x3),y3, title="Fresnel")
p4 = Plots.plot(x,y4, label=L"\mathbf{H}_1(x)", title="Struve")

Plots.plot(p1,p2,p3,p4,layout=(2,2))

Plots.savefig("./FewSpecialFunctionsTest/combinedplot.pdf")
