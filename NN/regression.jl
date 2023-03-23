using CSV, GLM, Plots, TypedTables

# use CSV package to import data from CSV file

data = CSV.File("NN/housingdata.csv")

X = data.size

Y = round.(Int, data.price / 1000)

t = Table(X = X, Y = Y)

# use Plots package to generate scatter plot of data

gr(size = (600, 600))

# create scatter plot

p_scatter = scatter(X, Y,
    xlims = (0, 5000),
    ylims = (0, 800),
    xlabel = "Size (sqft)",
    ylabel = "Price (in thousands of dollars)",
    title = "Housing Prices in Portland",
    legend = false,
    color = :red
)

# use GLM package for Linear Regression model

ols = lm(@formula(Y ~ X), t)

# add linear regression line to plot

plot!(X, predict(ols), color = :green, linewidth = 3)

# predict price based on a new value for size

newX = Table(X = [1250])

predict(ols, newX)