using Flux, Plots, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10
)
scalefontsizes()

x = hcat(collect(Float32, -3:0.1:3)...)

f(x) = @. 3x + 2;

y = f(x)

x = x .* reshape(rand(Float32, 61), (1, 61))

plot(vec(x), vec(y), lw = 3, seriestype = :scatter, label = "", title = "Generated data", xlabel = "x", ylabel= "y")


#Build model
custom_model(W, b, x) = @. W*x + b;

W = rand(Float32, 1, 1)
b = [0.0f0]
custom_model(W, b, x) |> size
custom_model(W, b, x)[1], y[1]

function custom_loss(W, b, x, y)
    ŷ = custom_model(W, b, x)
    sum((y .- ŷ).^2) / length(x)
end;
custom_loss(W, b, x, y)
flux_model = Dense(1 => 1)
flux_model.weight, flux_model.bias

flux_model(x)[1], y[1]
function flux_loss(flux_model, x, y)
    ŷ = flux_model(x)
    Flux.mse(ŷ, y)
end;
flux_loss(flux_model, x, y)
W = Float32[1.1412252]
custom_loss(W, b, x, y), flux_loss(flux_model, x, y)

dLdW, dLdb, _, _ = gradient(custom_loss, W, b, x, y);
W .= W .- 0.1 .* dLdW

b .= b .- 0.1 .* dLdb
custom_loss(W, b, x, y)

function train_custom_model()
    dLdW, dLdb, _, _ = gradient(custom_loss, W, b, x, y)
    @. W = W - 0.1 * dLdW
    @. b = b - 0.1 * dLdb
end;
train_custom_model();

W, b, custom_loss(W, b, x, y)

for i = 1:40
    train_custom_model()
 end

 W, b, custom_loss(W, b, x, y)

 plot(reshape(x, (61, 1)), reshape(y, (61, 1)), lw = 3, seriestype = :scatter, label = "Data", title = "Simple Linear Regression", xlabel = L"x", ylabel= "Testing the font");
 plot!((x) -> b[1] + W[1] * x, -3, 3, label="Custom model", lw=2)  

 savefig("./figure7.pdf")
