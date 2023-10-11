# Import required libraries
using ITensors
using ITensorNetworks

# Define the number of spins and the bond dimension
num_spins = 10
bond_dimension = 5

# Define the single-site operators for the spin chain
Sz = diagonal([1.0, -1.0])
Sx = 0.5 * (Sz + diagm([1.0, -1.0]))

# Define the Hamiltonian of the spin chain (Ising model)
hamiltonian = sum(Sz(i) * Sz(i + 1) + 0.5 * Sx(i) for i in 1:num_spins - 1)

# Create random MPS tensors for the initial state
init_state = [randomITensor(bond(i-1), spin(i), bond(i)) for i in 1:num_spins]

# Define the MPS network
mps_network = MPS(init_state)

# Define the energy function for the variational method
function energy_fn(network::MPS)
    # Apply the Hamiltonian to the MPS network
    energy = overlap(network, hamiltonian, network)
    
    return energy
end

# Perform variational optimization to find the ground state
ground_state, energy = dmrg(mps_network, hamiltonian)

# Print the energy of the ground state
println("Ground state energy: ", energy)
