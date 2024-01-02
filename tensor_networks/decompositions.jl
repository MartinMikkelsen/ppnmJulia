using LinearAlgebra, TensorOperations

d1 = 10; d2 = 6;

A = rand(d1, d2);
F = svd(A);

Af = F.U*diagm(0 => F.S)*F.Vt;
dA 