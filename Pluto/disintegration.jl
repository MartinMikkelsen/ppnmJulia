using LinearAlgebra, Optim, Plots

RND() = rand() * 2 - 1
NORM(a) = 1

function print_matrix(A)
    println(A)
end

function print_vector(v)
    println(v)
end

function main()
    n1 = 1
    n2 = 5
    nt = n1 + n2
    hbarc = 197.3
    mpi = 139
    mN = 939
    mr = mpi * mN / (mpi + mN)
    bw = 2
    bscale = 3 * bw
    kappa = 1 / bw / bw
    A = 10 / bw
    H = zeros(nt, nt)
    N = zeros(nt, nt)

    function master(a)
        H[1, 1] = 0
        N[1, 1] = 1
        for k in 1:n2
            i = k + n1
            alpha_k = a[k]
            h0i = 3 * A * NORM(alpha_k) * 1.5 / (alpha_k + kappa) * (π / (alpha_k + kappa)) ^ 1.5
            H[1, i] = h0i
            H[i, 1] = h0i
            N[1, i] = 0
            N[i, 1] = 0
        end
        for k in 1:n2
            for l in k:n2
                i = k + n1
                j = l + n1
                alpha_k = a[k]
                alpha_l = a[l]
                normij = 3 * NORM(alpha_k) * NORM(alpha_l) * 1.5 / (alpha_k + alpha_l) * (π / (alpha_k + alpha_l)) ^ 1.5
                hamij = 3 * NORM(alpha_k) * NORM(alpha_l) * hbarc * hbarc / 2 / mr * 15 * alpha_k * alpha_l / (alpha_k + alpha_l) ^ 2 * (π / (alpha_k + alpha_l)) ^ 1.5 + mpi * normij
                H[i, j] = hamij
                H[j, i] = hamij
                N[i, j] = normij
                N[j, i] = normij
            end
        end
        F = eigen(H, N)
        E = F.values
        V = F.vectors
        sort_indices = sortperm(E)
        E = E[sort_indices]
        V = V[:, sort_indices]
        E0 = E[1]
        return E0, V
    end

    start = zeros(n2)
    step = ones(n2)
    for i in 1:n2
        b = bscale * (i + 1) / n2
        start[i] = 1 / b / b
    end

    E0, V = master(start)
    println("E0 = ", E0)

    result = optimize(master, start, NelderMead())
    println(result)

    function phi(r, result_x, V)
        sum = 0
        sign = sign(V[n1 + 1, 1])
        for i in 1:n2
            ci = V[i + n1, 1]
            alpha_i = result_x[i]
            sum += ci * exp(-alpha_i * r * r) * NORM(alpha_i)
        end
        return sign * sum
    end

    rmax = 6
    r = 0:1/16:rmax
    phi_r = phi.(r, result.minimizer, V)
    println("phi = ", phi_r)
    for (r_val, phi_val) in zip(r, phi_r)
        println(r_val, " ", phi_val)
    end
    plot(r, phi_r)
end

main()
