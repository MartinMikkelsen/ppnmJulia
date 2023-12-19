using DynamicalSystems, Plots


plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10
)


hh = Systems.henonheiles()
# create some initial conditions, all at energy = 0.13
ics = Systems.henonheiles_ics(0.13, 15)
for ic in ics
    psos = poincaresos(hh, (1, 0.0), 2000; u0 = ic) 
    λ = lyapunov(hh, 10000; u0 = ic)
    v = clamp(λ/0.06, 0, 1)
    p = plot!(psos[:, 2], psos[:, 4])
    display(p)
end
