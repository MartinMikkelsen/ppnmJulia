using LinearAlgebra, TensorOperations

d = 12;
A = zeros(d,d,d);
for ni = 1:d
    for nj = 1:d
        for nk = 1:d
            A[ni,nj,nk] = sqrt(ni + 2*nj + 3*nk);
        end
    end
end
B = A;
C = A;
H = ncon(Any[A,B,C],Any[[-1,-2,1],[1,-3,2],[2,-4,-5]]);
nH = norm(H)

um,sm,vm = svd(reshape(C,d,d^2));
chi = 2;
CL = um[:,1:chi]*diagm(0 => sqrt.(sm[1:chi]));
CR = reshape(diagm(0 => sqrt.(sm[1:chi]))*(vm[:,1:chi]'),chi,d,d);
H1 = ncon(Any[A,B,CL,CR],Any[[-1,-2,1],[1,-3,2],[2,3],[3,-4,-5]]);
ϵ = norm(H1-H)/nH

# iterate QR decomps
Qm, Rm = qr(reshape(A,d^2,d));
Ap = reshape(Matrix(Qm),d,d,d);
Qm, Rm = qr(reshape(ncon(Any[Rm,B],Any[[-1,1],[1,-2,-3]]),d^2,d));
Bp = reshape(Matrix(Qm),d,d,d);
Cp = ncon(Any[Rm,C],Any[[-1,1],[1,-2,-3]]);
# check result
Hp = ncon(Any[Ap,Bp,Cp],Any[[-1,-2,1],[1,-3,2],[2,-4,-5]]);
errp = norm(H - Hp) / nH; # should be zero

# (d) globally optimal truncated SVD
um,sm,vm = svd(reshape(Cp,d,d^2));
chi = 2;
CLp = um[:,1:chi]*diagm(0 => sqrt.(sm[1:chi]));
CRp = reshape(diagm(0 => sqrt.(sm[1:chi]))*(vm[:,1:chi]'),chi,d,d);

H1p = ncon(Any[Ap,Bp,CLp,CRp],Any[[-1,-2,1],[1,-3,2],[2,3],[3,-4,-5]]);
err1p = norm(H - H1p) / nH;

