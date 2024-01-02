using LinearAlgebra, TensorOperations

d = 5; # local dimension
chi = 3; # max internal dimension
H0 = reshape(sqrt.(1:d^7),d,d,d,d,d,d,d); # initial tensor

# first decomposition
utemp,stemp,vtemp = svd(reshape(H0,d^2,d^5));
U0 = reshape(utemp[:,1:chi],d,d,chi);
H1 = reshape(diagm(0 => stemp[1:chi])*vtemp[:,1:chi]',chi,d,d,d,d,d);
# second decomposition
utemp,stemp,vtemp = svd(reshape(permutedims(H1,[2,3,1,4,5,6]),d^2,chi*d^3));
U1 = reshape(utemp[:,1:chi],d,d,chi);
H2 = permutedims(reshape(diagm(0 => stemp[1:chi])*vtemp[:,1:chi]',chi,chi,d,d,d),[2,1,3,4,5]);
# third decomposition
utemp,stemp,vtemp = svd(reshape(H2,chi^2,d^3));
U2 = reshape(utemp[:,1:chi],chi,chi,chi);
H3 = reshape(diagm(0 => stemp[1:chi])*vtemp[:,1:chi]',chi,d,d,d);
# fourth decomposition
utemp,stemp,vtemp = svd(reshape(H3,chi*d,d^2));
V3 = reshape(conj(vtemp[:,1:chi]),d,d,chi);
H4 = reshape(utemp[:,1:chi]*diagm(0 => stemp[1:chi]),chi,d,chi);
# check result
H0recovered = ncon(Any[U0,U1,U2,V3,H4],Any[[-1,-2,1],[-3,-4,2],[1,2,3],[-6,-7,4],[3,-5,4]]);
totErr = norm(H0 - H0recovered) / norm(H0)

##### Ex4.2(a): set B-C link as center of orthogonality
d = 5; # index dimension
A = rand(d,d,d);
B = rand(d,d,d);
C = rand(d,d,d);
Sig = Matrix{Float64}(I, d, d); # initial link matrix

# generate gauge change matrices
rho1 = ncon(Any[A,A,B,B],Any[[1,2,3],[1,2,4],[3,5,-1],[4,5,-2]]);
rho2 = ncon(Any[C,C],Any[[-1,1,2],[-2,1,2]]);
d1, u1 = eigen(rho1); sq_d1 = sqrt.(abs.(d1));
d2, u2 = eigen(rho2); sq_d2 = sqrt.(abs.(d2));
X1 = u1*diagm(0 => sq_d1)*u1'; X1inv = u1*diagm(0 => 1 ./sq_d1)*u1';
X2 = u2*diagm(0 => sq_d2)*u2'; X2inv = u2*diagm(0 => 1 ./sq_d2)*u2';
# implement gauge change
Bprime = ncon(Any[B,X1inv],Any[[-1,-2,1],[1,-3]]);
Cprime = ncon(Any[X2inv,C],Any[[-1,1],[1,-2,-3]]);
Sig_prime = X1*Sig*X2;
# check result
H0 = ncon(Any[A,B,C],Any[[-1,-2,1],[1,-3,2],[2,-4,-5]]);
H1 = ncon(Any[A,Bprime,Sig_prime,Cprime],Any[[-1,-2,1],[1,-3,2],[2,3],[3,-4,-5]]);
totErr = norm(H0 - H1) / norm(H0)

# perform unitary gauge change to diagonalize link matrix
utemp, Sig_pp, vtemp = svd(Sig_prime);
Bpp = ncon(Any[Bprime,utemp],Any[[-1,-2,1],[1,-3]]);
Cpp = ncon(Any[Cprime,vtemp],Any[[1,-2,-3],[1,-1]]);
# check result
H2 = ncon(Any[A,Bpp,diagm(0 => Sig_pp),Cpp],Any[[-1,-2,1],[1,-3,2],[2,3],[3,-4,-5]]);
totErr = norm(H0 - H2) / norm(H0)

H = ncon(Any[A,B,C],Any[[-1,-2,1],[1,-3,2],[2,-4,-5]]);
utemp,stemp,vtemp = svd(reshape(H,d^3,d^2));
UH = reshape(utemp[:,1:d],d,d,d,d);
SH = diagm(0 => stemp[1:d]);
VH = reshape(vtemp[:,1:d]',d,d,d);

# compare with previous tensors from orthonormal form
ErrU = norm(abs.(UH[:]) - abs.(reshape(ncon(Any[A,Bpp],Any[[-1,-2,1],[1,-3,-4]]),d^4)))
ErrS = norm(diag(SH) - Sig_pp)
ErrV = norm(abs.(VH[:]) - abs.(Cpp[:]))
# all three results should be vanishingly small!!!