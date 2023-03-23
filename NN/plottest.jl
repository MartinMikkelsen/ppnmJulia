using Flux, Plots
data = [([x], 13x-x^7) for x in -2:0.1f0:2]

model = Chain(Dense(1 => 23, tanh), Dense(23 => 1, bias=false), only)

optim = Flux.setup(Adam(), model)
for epoch in 1:5000
  Flux.train!((m,x,y) -> (m(x) - y)^2, model, data, optim)
end

plot(x -> 13x-x^7, -2, 2, legend=false)
scatter!(-2:0.1:2, [model([x]) for x in -2:0.1:2])