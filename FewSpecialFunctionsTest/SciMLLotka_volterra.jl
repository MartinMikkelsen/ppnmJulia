using DifferentialEquations, Optimization, OptimizationPolyalgorithms, SciMLSensitivity
using ForwardDiff, Plots

function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α * x - β * x * y
    du[2] = dy = -δ * y + γ * x * y
end

# Initial condition
u0 = [1.0, 1.0]

# Simulation interval
tspan = (0.0, 10.0)

# LV equation parameter. p = [α, β, δ, γ]
p = [1.5, 1.0, 3.0, 1.0]

# Setup the ODE problem, then solve
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
datasol = solve(prob, saveat = 1)
data = Array(datasol)

## Now do the optimization process
function loss(newp)
    newprob = remake(prob, p = newp)
    sol = solve(newprob, saveat = 1)
    loss = sum(abs2, sol .- data)
    return loss, sol
end

callback = function (p, l, sol)
    display(l)
    plt = plot(sol, ylim = (0, 6), label = "Current Prediction")
    scatter!(plt, datasol, label = "Data")
    display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = Optimization.AutoForwardDiff()
pguess = [1.0, 1.2, 2.5, 1.2]
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, pguess)

result_ode = Optimization.solve(optprob, PolyOpt(),
                                callback = callback,
                                maxiters = 200)
