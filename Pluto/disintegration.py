import numpy as np
from scipy import linalg, optimize
import matplotlib.pyplot as plt

RND = lambda: np.random.rand() * 2 - 1
NORM = lambda a: 1

def print_matrix(A):
    print(A)

def print_vector(v):
    print(v)

def main():
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
    H = np.zeros((nt, nt))
    N = np.zeros((nt, nt))
    copyN = N.copy()
    copyH = H.copy()

    def master(a, params=None):
        nonlocal H, N
        H[0, 0] = 0
        N[0, 0] = 1
        for k in range(n2):
            i = k + n1
            alpha_k = a[k]
            h0i = 3 * A * NORM(alpha_k) * 1.5 / (alpha_k + kappa) * pow(np.pi / (alpha_k + kappa), 1.5)
            H[0, i] = h0i
            H[i, 0] = h0i
            N[0, i] = 0
            N[i, 0] = 0
        for k in range(n2):
            for l in range(k, n2):
                i = k + n1
                j = l + n1
                alpha_k = a[k]
                alpha_l = a[l]
                normij = 3 * NORM(alpha_k) * NORM(alpha_l) * 1.5 / (alpha_k + alpha_l) * pow(np.pi / (alpha_k + alpha_l), 1.5)
                hamij = 3 * NORM(alpha_k) * NORM(alpha_l) * hbarc * hbarc / 2 / mr * 15 * alpha_k * alpha_l / pow(alpha_k + alpha_l, 2) * pow(np.pi / (alpha_k + alpha_l), 1.5) + mpi * normij
                H[i, j] = hamij
                H[j, i] = hamij
                N[i, j] = normij
                N[j, i] = normij
        E = linalg.eigh(H, N, eigvals_only=False)[0]
        V = linalg.eigh(H, N, eigvals_only=False)[1]
        idx = E.argsort()
        E = E[idx]
        V = V[:, idx]
        E0 = E[0]
        return E0

    start = np.zeros(n2)
    step = np.zeros(n2)
    for i in range(n2):
        b = bscale * (i + 1) / n2
        start[i] = 1 / b / b
        step[i] = 1

    E0 = master(start)
    print(f"E0={E0}")

    result = optimize.minimize(master, start, method="Nelder-Mead")
    print(result)

    def phi(r):
        sum = 0
        V = linalg.eigh(H, N, eigvals_only=False)[1]
        sign = np.sign(V[n1, 0])
        for i in range(n2):
            ci = V[i + n1, 0]
            alpha_i = result.x[i]
            sum += ci * np.exp(-alpha_i * r * r) * NORM(alpha_i)
        return sign * sum

    rmax = 6
    r = np.arange(0, rmax, 1 / 16)
    phi_r = phi(r)
    print("phi=", phi_r)
    for r_val, phi_val in zip(r, phi_r):
        print(f"{r_val} {phi_val}")
    plt.plot(r,phi_r)
if __name__ == "__main__":
    main()
    plt.show()
