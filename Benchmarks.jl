using ParameterizedFunctions, OrdinaryDiffEq,
      ODEInterfaceDiffEq, Plots, Sundials, SciPyDiffEq, deSolveDiffEq
using DiffEqDevTools
using LinearAlgebra, StaticArrays

f = @ode_def_bare LotkaVolterra begin
    dx = a*x - b*x*y
    dy = -c*y + d*x*y
  end a b c d
  p = [1.5,1,3,1]
  tspan = (0.0,10.0)
  u0 = [1.0,1.0]
  prob = ODEProblem(f,u0,tspan,p)
  staticprob = ODEProblem{false}(f,SVector{2}(u0),tspan,SVector{4}(p))
  
  sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
  test_sol = TestSolution(sol)
  
  setups = [
            Dict(:alg=>DP5())
            Dict(:alg=>Tsit5())
            Dict(:alg=>Vern7())
            Dict(:prob_choice => 2, :alg=>DP5())
            Dict(:prob_choice => 2, :alg=>Tsit5())
            Dict(:prob_choice => 2, :alg=>Vern7())
            Dict(:alg=>dopri5())
            Dict(:alg=>SciPyDiffEq.RK45())
            Dict(:alg=>SciPyDiffEq.LSODA())
            Dict(:alg=>SciPyDiffEq.odeint())
            Dict(:alg=>deSolveDiffEq.lsoda())
            Dict(:alg=>deSolveDiffEq.ode45())
            Dict(:alg=>CVODE_Adams())
    ]
  
  labels = [
    "Julia: DP5"
    "Julia: Tsit5"
    "Julia: Vern7"
    "Julia: DP5 Static"
    "Julia: Tsit5 Static"
    "Julia: Vern7 Static"
    "Hairer: dopri5"
    "MATLAB: ode45"
    "MATLAB: ode113"
    "SciPy: RK45"
    "SciPy: LSODA"
    "SciPy: odeint"
    "deSolve: lsoda"
    "deSolve: ode45"
    "Sundials: Adams"
    ]
  
  abstols = 1.0 ./ 10.0 .^ (6:13)
  reltols = 1.0 ./ 10.0 .^ (3:10)
  wp = WorkPrecisionSet([prob,staticprob],abstols,reltols,setups;
                        names = labels,print_names = true,
                        appxsol=[test_sol,test_sol],dense=false,
                        save_everystep=false,numruns=100,maxiters=10000000,
                        timeseries_errors=false,verbose=false)
  plot(wp,title="Non-stiff 1: Lotka-Volterra",legend=:outertopleft,
       color=permutedims([repeat([:LightGreen],3)...,repeat([:DarkGreen],3)...,
       :Red,repeat([:Orange],2)...,repeat([:Yellow],3)...,
       repeat([:Blue],2)...,:Purple]),size = (800,350),
       xticks = 10.0 .^ (-12:1:5),
       yticks = 10.0 .^ (-6:0.5:5),
       bottom_margin=5Plots.mm)
  