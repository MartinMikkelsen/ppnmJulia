using LinearAlgebra, TensorOperations

d1 = 10; d2 = 8;

A = rand(d1,d1,d2,d2);

F = svd(reshape(A,d1^2,d2^2));
F.S

sqrt.(ncon(Any[A,conj(A)],Any[collect(1:ndims(A)),collect(1:ndims(A))]))

Aprime = A./norm(A);

Fprime = svd(reshape(Aprime,d1^2,d2^2));

sqrt(sum((Fprime.S).^2))

deltaval = 1e-4;
r_delta = sum(Fprime.S .> deltaval)
eps_err = sqrt.(sum(Fprime.S[r_delta:end] .^2))