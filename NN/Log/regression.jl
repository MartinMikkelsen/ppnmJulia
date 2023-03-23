using Plots, CSV

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10
)


theta_0 = 1.0
theta_1 = 1.0

z(x) = theta_0 .+ theta_1 *x

h(x) = 1 ./ (1+exp(-z(x)))

x = range(-10,10,100)

plot(x,h)

X = [0.247,0.327,0.938,0.347,0.299,0.245,0.285,0.360,0.394,0.509,1.045,0.638,0.816,1.036,0.594,0.398,0.363,0.561,0.432,0.409,0.569,0.473,0.529,0.356,0.421, 0.656,0.400,0.853]
Y_temp = ["absent","present","present","present","present","absent","present","present","present","present","present","present","present","present","present","absent","absent","absent","absent","absent","absent","present","present","absent","present","present","present","present"]

Y = []

for i in 1:length(Y_temp)
    if Y_temp[i] == "present"
        y = 1.0
    else 
        y = 0.0
    end
    push!(Y,y)
end

p_data = scatter(X,Y)

#initialize parameters

theta_0 = 0.0

theta_1 = 1.0

t0_history = []
t1_history = []

push!(t0_history, theta_0)
push!(t1_history, theta_1)

# Define hypothesis function
z(x) = theta_0 .+ theta_1*x
h(x) = 1 ./ (1 .+ exp.(-z(x))) 


plot(x,z)

# Define Cost function
m = length(X)

y_hat = h(X)

function cost()
    (-1 / m) * sum(Y .* log.(y_hat)) + (1 .- Y) .* log(1 .- y_hat)
end


