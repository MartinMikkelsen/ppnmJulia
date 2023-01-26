using DifferentialEquations

S = 1.0
b = 2.0
mu = 3.0
m_pi = 4.0

# Defining the differential equation
f(r,ϕ,dϕdr,E) = - (1 / (2 * mu)) * (dϕdr + (2 / r) * ϕ) + m_pi * ϕ - E * ϕ

# Defining the initial condition and boundary condition
initial_condition = S * exp(-r^2 / b^2)

# Defining the domain and the neural network architecture
u0 = initial_condition
tspan = (0.0, Inf)
p = [mu,m_pi]
prob = PDEProblem(f,u0,tspan,p)

# Solving the PDE
sol = solve(prob,Tsit5())

# Integrating the solution
E = 12 * pi * quadgk(x -> x^4 * sol(x)[1], 0, Inf)[1]

sol, E



