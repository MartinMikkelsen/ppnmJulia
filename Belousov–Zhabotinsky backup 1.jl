using DifferentialEquations, Plots, Term

s = 77.27
w = 0.161
q = 8.375*1e-6

function Bel_Zha(du,u,p,t)
    du[1] = s*(u[2]-u[1]*u[2]+u[1]-q*u[1]^2)
    du[2] = (-u[2]-u[1]*u[2]+u[3])./s
    du[3] = w*(u[1]-u[3])
end

u0 = [1.0, 2.0, 3.0]
tspan = (0.0,360.0)
p = [77.27,8.375e-6,0.161]
prob = ODEProblem(Bel_Zha,[1.0,2.0,3.0],(0.0,360.0),p)
sol = solve(prob)
plot(sol,linewidth=2.5,legend=:top)
xlabel!("t")
ylabel!("vals")
title!("Solutions")