using DynamicalSystems # load the library
#testing DynamicalSystems
function lorenz_rule(u, p, t) # the dynamics as a function
    σ, ρ, β = p
    x, y, z = u
    dx = σ*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - β*z
    return SVector(dx, dy, dz) # Static Vector
end
p  = [10.0, 28.0, 8/3] # parameters: σ, ρ, β
u0 = [0.0, 10.0, 0.0]  # initial condition
# create an instance of a `DynamicalSystem`
lorenz = ContinuousDynamicalSystem(lorenz_rule, u0, p)
T  = 100.0 # total time
Δt = 0.01  # sampling time
A  = trajectory(lorenz, T; Δt)