using LinearAlgebra, TensorOperations

d = 6;
H = zeros(d,d,d,d,d);
for ni = 1:d
    for nj = 1:d
        for nk = 1:d
            for nl = 1:d
                for nm = 1:d
                    H[ni,nj,nk,nl,nm] = sqrt(ni+2*nj+3*nk+4*nl+5*nm);
                end
            end
        end
    end
end
nmH = norm(H)
H1 = H/nmH;

# (b) multi-stage decomposition
chi = 4;
utemp,stemp,vtemp = svd(reshape(H1,d^2,d^3));
A = reshape(utemp[:,1:chi],d,d,chi);
Htemp = reshape(diagm(0 => stemp[1:chi])*(vtemp[:,1:chi]'),chi,d,d,d);
utemp,stemp,vtemp = svd(reshape(Htemp,chi*d,d^2));
B = reshape(utemp[:,1:chi],chi,d,chi);
C = reshape(diagm(0 => stemp[1:chi])*(vtemp[:,1:chi]'),chi,d,d);
# check accuracy
H2 = ncon(Any[A,B,C],Any[[-1,-2,1],[1,-3,2],[2,-4,-5]]);
TrErr = norm(H1-H2)

# (c) center of orthogonality at a link
# compute gauge change between A-B tensors
SigAB = Matrix{Float64}(I, chi, chi);
rhoL = ncon(Any[A,A],Any[[1,2,-1],[1,2,-2]]);
rhoR = ncon(Any[B,B,C,C],Any[[-1,1,2],[-2,1,3],[2,4,5],[3,4,5]]);
dL,uL = eigen(rhoL); sq_dL = sqrt.(abs.(dL));
dR,uR = eigen(rhoR); sq_dR = sqrt.(abs.(dR));
XL = uL*diagm(0 => sq_dL)*uL'; XLinv = uL*diagm(0 => 1 ./sq_dL)*uL';
XR = uR*diagm(0 => sq_dR)*uR'; XRinv = uR*diagm(0 => 1 ./sq_dR)*uR';
utemp,SigABp,vtemp = svd(XL*SigAB*XR);
Ap = ncon(Any[A,XLinv*utemp],Any[[-1,-2,1],[1,-3]]);
Bp = ncon(Any[B,vtemp'*XRinv],Any[[1,-2,-3],[-1,1]]);

# difference in singular values between SVD:
~,svalues,~ = svd(reshape(H2,d^2,d^3));
Serr = norm(svalues[1:chi] - SigABp)