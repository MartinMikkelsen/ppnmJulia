using LinearAlgebra, TensorOperations

d1 = 10; d2 = 6;

A = rand(d1, d2);
F = svd(A);

Af = F.U*diagm(0 => F.S)*F.Vt;
dA = norm(Af[:]-A[:])

### Higher order

d = 10; A = rand(d,d,d);
F = svd(reshape(A,d^2,d));
U = reshape(F.U,d,d,d);
# check result
Af = ncon(Any[U,diagm(0 => F.S),F.Vt],
    Any[[-1,-2,1],[1,2],[2,-3]]);
dA = norm(Af[:]-A[:])

##### Initialize unitaries and isometries
d1 = 10; d2 = 6;

# d1-by-d1 random unitary matrix U
U = svd(rand(d1,d1)).U

# d1-by-d2 random isometric matrix W
W = svd(rand(d1,d2)).U

##### Ex2.2(c): spect. decomp. of matrix
d = 10; A = rand(d,d);
H = 0.5*(A+A'); #random Hermitian
F = eigen(H);
# check result
Hf = F.vectors*diagm(0 => F.values)*F.vectors'
dH = norm(Hf[:]-H[:])

##### Ex2.2(d): spect. decomp. of tensor
d = 2; A = rand(d,d,d,d);
H = 0.5*(A + permutedims(A,[3,4,1,2]));
F = eigen(reshape(H,d^2,d^2));
U = reshape(F.vectors,d,d,d^2);
D = diagm(0 => F.values);
# check result
Hf = ncon(Any[U,D,U],
    Any[[-1,-2,1],[1,2],[-3,-4,2]]);
dH = norm(Hf[:]-H[:])

##### Ex2.2(f): QR decomp of matrix
d1 = 10; d2 = 6;
A = rand(d1,d2);
F = qr(A);
# check result
Af = F.Q*F.R;
dA = norm(Af[:]-A[:])

##### Ex2.2(g): QR decomp of tensor
d = 10;
A = rand(d,d,d);
F = qr(reshape(A,d^2,d));
Q = reshape(Array(F.Q),d,d,d);
# check result
Af = ncon(Any[Q,F.R],Any[[-1,-2,1],[1,-3]]);
dA = norm(Af[:]-A[:])

##### Ex2.3(c)
d = 10;
A = rand(10,10,10,10,10);
# frobenus norm
frobA0 = sqrt.(ncon(Any[A,conj(A)],Any[collect(1:ndims(A)),collect(1:ndims(A))]))
# equivalent frobenus norm
frobA1 = sqrt(sum(abs.(A[:]).^2))
# also equivalent frobenus norm
frobA2 = norm(A)