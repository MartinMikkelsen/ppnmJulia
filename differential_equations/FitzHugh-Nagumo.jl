using DynamicalSystems, OrdinaryDiffEq, Plots

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10
)

ds = Systems.fitzhugh_nagumo([0.0,0.0]; I = 0.0)
pulses_start = [20, 80, 170]
pulses_end = pulses_start .+ 4 # 4 = pulse width
pulses = sort!(vcat(pulses_start, pulses_end))
I = 0.2 # strength of pulses of I current
# Create the "callbacks": events in ODE solution 
condition(u,t,integ) = t ∈ pulses # trigger condition 
function affect!(integ) # what happens at the integrator
    i = integ.t ∈ pulses_start ? I : 0.0
    integ.p[4] = i # 4th parameter is value of current I
end
cb = DiscreteCallback(condition, affect!)
# transform `ds` to form allowing callbacks and solve:
prob = ODEProblem(ds, (0.0, 250.0))
sol = solve(prob, Tsit5(); callback=cb, tstops = pulses) 
plot(sol.t, sol[1, :]) # plot timeseries of u
pulse_ts = [any(x -> x ≤ t ≤ x+4, pulses_start) ? I : 0.0
            for t in sol.t]
plot(sol.t, pulse_ts, xlabel=L"t", ylabel=L"u")

