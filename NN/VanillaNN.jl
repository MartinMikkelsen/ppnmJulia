"""
The goal is to build a neural network without any packages in Julia
"""
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10
)

"""Make a function to initialize the model parameters. The parameters are the weight W and bias b"""

function initialise_model_weights(layer_dim, seed)
    params = Dict()

    for l=2:length(layer_dim)
        params[string("W_", (l-1))] = rand(StableRNG(seed), layer_dim[l], layer_dim[l-1]) * sqrt(2 / layer_dim[l-1])
        params[string("b_", (l-1))] = zeros(layer_dim[l], 1)
    end
    return params
end

"""Make a function that for some data makes a prediction through an activation function"""

function sigmoid(Z)
    A = 1 ./ (1 .+ exp.(.-Z))
    return (A = A, Z = Z)
end

function relu(Z)
    A = max.(0, Z)
    return (A = A, Z = Z)
end

function linear_forward(A,W,b)
    Z = (W*A).+b
    cache = (A,W,b)

    @assert size(Z) == (size(W,1),size(A,2))

    return (Z=Z, cache=cache)
end

function linear_forward_activation(A_prev, W, b, activation_function="relu")
    @assert activation_function ∈ ("sigmoid", "relu")
    Z, linear_cache = linear_forward(A_prev, W, b)

    if activation_function == "sigmoid"
        A, activation_cache = sigmoid(Z)
    end

    if activation_function == "relu"
        A, activation_cache = relu(Z)
    end

    cache = (linear_step_cache=linear_cache, activation_step_cache=activation_cache)

    @assert size(A) == (size(W, 1), size(A_prev, 2))

    return A, cache
end

function forward_propagate_model_weights(DMatrix, parameters)
    master_cache = []
    A = DMatrix
    L = Int(length(parameters) / 2)

    # Forward propagate until the last (output) layer
    for l = 1 : (L-1)
        A_prev = A
        A, cache = linear_forward_activation(A_prev,
                                             parameters[string("W_", (l))],
                                             parameters[string("b_", (l))],
                                             "relu")
        push!(master_cache , cache)
    end

    # Make predictions in the output layer
    Ŷ, cache = linear_forward_activation(A,
                                         parameters[string("W_", (L))],
                                         parameters[string("b_", (L))],
                                         "sigmoid")
    push!(master_cache, cache)

    return Ŷ, master_cache
end

function calculate_cost(Ŷ, Y)
    m = size(Y, 2)
    epsilon = eps(1.0)

    # Deal with log(0) scenarios
    Ŷ_new = [max(i, epsilon) for i in Ŷ]
    Ŷ_new = [min(i, 1-epsilon) for i in Ŷ_new]

    cost = -sum(Y .* log.(Ŷ_new) + (1 .- Y) .* log.(1 .- Ŷ_new)) / m
    return cost
end

function sigmoid_backwards(∂A,activated_cache)
    s = sigmoid(activated_cache).A
    ∂Z = ∂A .* s .*(1 .- s)

    @assert (size(∂Z) == size(activated_cache))
    return ∂Z
end

function relu_backwards(∂A,activated_cache)
    return ∂A .* (activated_cache .> 0)
end

function linear_backwards(∂Z, cache)
    A_prev, W, b = cache
    m = size(A_prev,2)

    ∂W = ∂Z .* (A_prev')/m
    ∂b = sum(∂Z, dim=2)/m
    ∂A_prev = (W')*∂Z

    @assert (size(∂A_prev) == size(A_prev))
    @assert (size(∂W) == size(W))
    @assert (size(∂b) == size(b))

    return ∂W, ∂b, ∂A_prev
end

function linear_backward(∂Z, cache)
    # Unpack cache
    A_prev , W , b = cache
    m = size(A_prev, 2)

    # Partial derivates of each of the components
    ∂W = ∂Z * (A_prev') / m
    ∂b = sum(∂Z, dims = 2) / m
    ∂A_prev = (W') * ∂Z

    @assert (size(∂A_prev) == size(A_prev))
    @assert (size(∂W) == size(W))
    @assert (size(∂b) == size(b))

    return ∂W , ∂b , ∂A_prev
end

"""
    Compute the gradients (∇) of the parameters (master_cache) of the constructed model
    with respect to the cost of predictions (Ŷ) in comparison with actual output (Y).
"""
function back_propagate_model_weights(Ŷ, Y, master_cache)
    # Initiate the dictionary to store the gradients for all the components in each layer
    ∇ = Dict()

    L = length(master_cache)
    Y = reshape(Y , size(Ŷ))

    # Partial derivative of the output layer
    ∂Ŷ = (-(Y ./ Ŷ) .+ ((1 .- Y) ./ ( 1 .- Ŷ)))
    current_cache = master_cache[L]

    # Backpropagate on the layer preceeding the output layer
    ∇[string("∂W_", (L))], ∇[string("∂b_", (L))], ∇[string("∂A_", (L-1))] = linear_activation_backward(∂Ŷ,
                                                                                                       current_cache,
                                                                                                       "sigmoid")
    # Go backwards in the layers and compute the partial derivates of each component.
    for l=reverse(0:L-2)
        current_cache = master_cache[l+1]
        ∇[string("∂W_", (l+1))], ∇[string("∂b_", (l+1))], ∇[string("∂A_", (l))] = linear_activation_backward(∇[string("∂A_", (l+1))],
                                                                                                             current_cache,
                                                                                                             "relu")
    end

    # Return the gradients of the network
    return ∇
end

"""
    Check the accuracy between predicted values (Ŷ) and the true values(Y).
"""
function assess_accuracy(Ŷ , Y)
    @assert size(Ŷ) == size(Y)
    return sum((Ŷ .> 0.5) .== Y) / length(Y)
end

"""
    Update the paramaters of the model using the gradients (∇)
    and the learning rate (η).
"""
function update_model_weights(parameters, ∇, η)

    L = Int(length(parameters) / 2)

    # Update the parameters (weights and biases) for all the layers
    for l = 0: (L-1)
        parameters[string("W_", (l + 1))] -= η .* ∇[string("∂W_", (l + 1))]
        parameters[string("b_", (l + 1))] -= η .* ∇[string("∂b_", (l + 1))]
    end

    return parameters
end

"""
    Train the network using the desired architecture that best possible
    matches the training inputs (DMatrix) and their corresponding ouptuts(Y)
    over some number of iterations (epochs) and a learning rate (η).
"""
function train_network(layer_dim , DMatrix, Y;  η=0.001, epochs=1000, seed=2020, verbose=true)
    # Initiate an empty container for cost, iterations, and accuracy at each iteration
    costs = []
    iters = []
    accuracy = []

    # Initialise random weights for the network
    params = initialise_model_weights(layer_dim, seed)

    # Train the network
    for i = 1:epochs

        Ŷ , caches  = forward_propagate_model_weights(DMatrix, params)
        cost = calculate_cost(Ŷ, Y)
        acc = assess_accuracy(Ŷ, Y)
        ∇  = back_propagate_model_weights(Ŷ, Y, caches)
        params = update_model_weights(params, ∇, η)

        if verbose
            println("Iteration -> $i, Cost -> $cost, Accuracy -> $acc")
        end

        # Update containers for cost, iterations, and accuracy at the current iteration (epoch)
        push!(iters , i)
        push!(costs , cost)
        push!(accuracy , acc)
    end
        return (cost = costs, iterations = iters, accuracy = accuracy, parameters = params)
end

using Plots
using MLJBase
using StableRNGs # Seeding generator for reproducibility

# Generate fake data
X, y = make_blobs(10_000, 3; centers=2, as_table=false, rng=2020);
X = Matrix(X');
y = reshape(y, (1, size(X, 2)));
f(x) =  x == 2 ? 0 : x
y2 = f.(y);

# Input dimensions
input_dim = size(X, 1);

# Train the model
nn_results = train_network([input_dim, 5, 3, 1], X, y2; η=0.01, epochs=50, seed=1, verbose=true);

# Plot accuracy per iteration
p1 = plot(nn_results.accuracy,
         label="Accuracy",
         xlabel="Number of iterations",
         ylabel="Accuracy as %",
         title="Development of accuracy at each iteration");

# Plot cost per iteration 
p2 = plot(nn_results.cost,
         label="Cost",
         xlabel="Number of iterations",
         ylabel="Cost (J)",
         color="red",
         title="Development of cost at each iteration");

# Combine accuracy and cost plots
plot(p1, p2, layout = (2, 1), size = (800, 600))

savefig("NN/Accuracy_Cost.pdf")
