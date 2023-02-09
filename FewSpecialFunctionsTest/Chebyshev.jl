using Plots, FewSpecialFunctions, SpecialFunctions, QuadGK

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)


function T_0(x)
    return 1
end
function T_1(x)
    return x
end
function T(n,x)
    return 2*x*T
end