using LinearAlgebra

function calculate_mk(masses::Vector{Float64})
    cumsum(masses)
end

function jacobi_matrix(masses::Vector{Float64})
    N = length(masses)
    M = calculate_mk(masses)
    Omega = zeros(N, N)
    for i in 1:N
        if i == 1
            Omega[i, i] = 1
            Omega[i, i+1] = -1
        elseif i < N
            for j in 1:i
                Omega[i, j] = masses[j]/M[i]
            end
            Omega[i, i+1] = -1
        else
            for j in 1:N
                Omega[i, j] = masses[j]/M[N]
            end
        end
    end
    Omega
end

function relative_coordinates_matrix(masses::Vector{Float64})
    N = length(masses)
    M = calculate_mk(masses)
    Omega = zeros(N, N)
    for i in 1:N
        if i < N
            Omega[i, i] = -1
            Omega[i, i+1] = 1
        else
            for j in 1:N
                Omega[i, j] = masses[j]/M[N]
            end
        end
    end
    Omega
end

function transform_coordinates(r::Matrix{Float64}, Omega::Matrix{Float64})
    Omega \ r
end

function transform_back(x::Matrix{Float64}, Omega::Matrix{Float64})
    Omega * x
end

# Usage
N = 5
particles = [rand(3) for _ in 1:N]
masses = rand(N)
jacobi_coords = relative_coordinates_matrix(masses)
Ω = jacobi_matrix(masses)
