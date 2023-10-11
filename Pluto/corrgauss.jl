using LinearAlgebra, Plots

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5,
    framestyle=:box,
    label=nothing,
    grid=true,
    palette=:tab10,
)


function correlated_gaussian(x, y, alpha, beta)
    return exp(-alpha*beta/(alpha+beta)*dot(x-y, x-y))
end

function overlap_matrix(a, b, alpha, beta, n)
    S = zeros(n,n)
    for i=1:n
        for j=1:n
            A = alpha[i] + alpha[j]
            B = beta[i] + beta[j]
            P = (alpha[i]*a[i,:] + alpha[j]*a[j,:])/A
            Q = (beta[i]*b[i,:] + beta[j]*b[j,:])/B
            S[i,j] = correlated_gaussian(P,Q,A,B,alpha[i],beta[i])*(pi/(A+B))^(3/2)
        end
    end
    return S
end

function potential_matrix(a, b, alpha, beta, n, V)
    V_mat = zeros(n,n)
    for i=1:n
        for j=1:n
            for k=1:n
                for l=1:n
                    A = alpha[i] + alpha[j]
                    B = beta[k] + beta[l]
                    P = (alpha[i]*a[i,:] + alpha[j]*a[j,:])/A
                    Q = (beta[k]*b[k,:] + beta[l]*b[l,:])/B
                    V_mat[i,j] += a[i,:]'*a[j,:]*V(P)*correlated_gaussian(P,Q,A,B,alpha[i],beta[k])*(pi/(A+B))^(3/2)
                end
            end
        end
    end
    return V_mat
end

function solve_eigenproblem(S, H)
    F = H
    eigvals, eigvecs = eigen(S)
    X = eigvecs
    C = X'*F*X
    epsilon, C_canonical = eigen(C)
    C = X*C_canonical
    return eigvals, C
end

a = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
b = [0.5 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.5]
alpha = [4, 3, 2.1]
beta = [0.9, 0.8, 0.3]
V(x) = 0.5*(x[1]-x[2])^2

n = length(alpha)
S = overlap_matrix(a, b, alpha, beta, n)
H = potential_matrix(a, b, alpha, beta, n, V)
eigvals, eigvecs = solve_eigenproblem(S, H)
phi = eigvecs[:, 1]
println("Eigenvalues: ", eigvals)

maxN = 25
eigvals_array = zeros(maxN)
for N=1:maxN
    alpha = [0.5]*N
    beta = [0.1]*N
    n = length(alpha)
    S = overlap_matrix(a, b, alpha, beta, n)
    H = potential_matrix(a, b, alpha, beta, n, V)
    eigvals, eigvecs = solve_eigenproblem(S, H)
    eigvals_array[N] = eigvals[1]
    phi = eigvecs[:, 1]
end

plot(1:maxN, eigvals_array, xlabel="Number of basis functions", ylabel="Ground state energy", legend=false)

plot(x -> begin
        y = [x, 0.0, 0.0]
        phi[1]*correlated_gaussian(a[1,:], y, alpha[1], 0.1, alpha[1], 0.1) + 
        phi[2]*correlated_gaussian(a[2,:], y, alpha[2], 0.1, alpha[2], 0.1) +
        phi[3]*correlated_gaussian(a[3,:], y, alpha[3], 0.1, alpha[3], 0.1)
    end,
    -5, 5, xlabel="x", ylabel="ϕ(x)")


phi_solution1(x) = (sqrt(2)/π)^(0.25)*exp(-1/(2*sqrt(2))*x^2)


plot!(x,phi_solution1.(x))