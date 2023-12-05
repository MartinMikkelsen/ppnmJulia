using LinearAlgebra
using Optim
using Plots

# Constants
n1, n2 = 1, 5
nt = n1 + n2
hbarc = 197.3
mpi = 139
mN = 939
mr = mpi * mN / (mpi + mN)
bw = 2
bscale = 3 * bw
kappa = 1 / bw / bw
A = 10 / bw

function calculate_elements(a)
    H = zeros(nt, nt)
    N = zeros(nt, nt)
    H[1, 1] = 0
    N[1, 1] = 1

    for k in 1:n2
        i = k + n1
        alpha_k = a[k]
        h0i = 3 * A * 1.5 / (alpha_k + kappa) * (π / (alpha_k + kappa)) ^ 1.5
        H[1, i] = h0i
        H[i, 1] = h0i
        N[1, i] = 0
        N[i, 1] = 0
    end

    for k in 1:n2, l in k:n2
        i, j = k + n1, l + n1
        alpha_k, alpha_l = abs(a[k]), abs(a[l])
        normij = 3 * 1.5 / (alpha_k + alpha_l) * (π / (alpha_k + alpha_l)) ^ 1.5
        hamij = (3 * hbarc^2 / 2 / mr * 15 * alpha_k * alpha_l / (alpha_k + alpha_l)^2 * (π / (alpha_k + alpha_l)) ^ 1.5) + mpi * normij
        H[i, j] = hamij
        H[j, i] = hamij
        N[i, j] = normij
        N[j, i] = normij
    end

    return H, N
end

# Objective function to minimize
function master(a)
    H, N = calculate_elements(a)
    E, V = eigen(Hermitian(H), Hermitian(N))
    E0 = minimum(E)
    return E0
end
start = [(bscale * (i + 1) / n2)^-2 for i in 1:n2]

results = optimize(master, start, NelderMead(), Optim.Options(iterations=100, x_tol=1e-5, f_tol=1e-5))
a = optimized_params = Optim.minimizer(results)
E0_optimized = Optim.minimum(results)

println("Optimized E0 = $E0_optimized")

_, V_optimized = calculate_elements(optimized_params)

# Define the phi function
function phi(r, V, a)
    sum = 0.0
    sign_V = sign(V[n1 + 1, 1])
    for i in 1:n2
        ci = V[i + n1, 1]
        alpha_i = a[i]
        sum += ci * exp(-alpha_i * r^2)
    end
    return sign_V * sum
end

rmax = 6
for r in 0:1/16:rmax
    println("$r $(phi(r, V_optimized, optimized_params))")
end

r_values = 0:0.0625:rmax  # Adjust step size as needed
phi_values = [phi(r, V, a) for r in r_values]

# Create the plot
plot(r_values, phi_values, label="phi(r)", xlabel="r", ylabel="phi", title="Wave Function")
