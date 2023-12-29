using LinearAlgebra, TensorOperations

##### Ex.3.3(c): Creating a center of orthogonality by 'pulling through'
# define tensors
d = 3;
A = rand(d,d,d,d); B = rand(d,d,d);
C = rand(d,d,d); D = rand(d,d,d);
E = rand(d,d,d); F = rand(d,d,d);
G = rand(d,d,d);
# iterate QR decomps
DQ, DR = qr(reshape(D,d^2,d)); DQ = reshape(Matrix(DQ),d,d,d);
EQ, ER = qr(reshape(E,d^2,d)); EQ = reshape(Matrix(EQ),d,d,d);
Btilda = ncon(Any[B,DR,ER],Any[[1,2,-3],[-1,1],[-2,2]]);
BQ, BR = qr(reshape(Btilda,d^2,d)); BQ = reshape(Matrix(BQ),d,d,d);
FQ, FR = qr(reshape(F,d^2,d)); FQ = reshape(Matrix(FQ),d,d,d);
GQ, GR = qr(reshape(G,d^2,d)); GQ = reshape(Matrix(GQ),d,d,d);
Ctilda = ncon(Any[C,GR],Any[[1,-2,-3],[-1,1]]);
CQ, CR = qr(reshape(Ctilda,d^2,d)); CQ = reshape(Matrix(CQ),d,d,d);
Aprime = ncon(Any[A,BR,FR,CR],Any[[1,-2,2,3],[-1,1],[-3,2],[-4,3]]);
# new network is formed from tensors: {Aprime,BQ,CQ,DQ,EQ,FQ,GQ}.

# check both networks evaluate to the same tensor
connectlist = Any[[3,-5,4,5],[1,2,3],[6,-10,5],[-1,-2,1],[-3,-4,2],[-6,-7,4],[-8,-9,6]];
H0 = ncon(Any[A,B,C,D,E,F,G],connectlist)[:];
H1 = ncon(Any[Aprime,BQ,CQ,DQ,EQ,FQ,GQ],connectlist)[:];
dH = norm(H0-H1)/norm(H0)

##### Ex.3.3(c): Creating a center of orthogonality with 'direct orthogonalization'
# define tensors
d = 3;
A = rand(d,d,d,d); B = rand(d,d,d);
C = rand(d,d,d); D = rand(d,d,d);
E = rand(d,d,d); F = rand(d,d,d);
G = rand(d,d,d);
# compute density matrices and their principle square roots
rho1 = ncon(Any[B,D,E,B,D,E],Any[[5,6,-2],[1,2,5],[3,4,6],[7,8,-1],[1,2,7],[3,4,8]]);
rho2 = ncon(Any[F,F],Any[[1,2,-2],[1,2,-1]]);
rho3 = ncon(Any[C,G,C,G],Any[[3,5,-2],[1,2,3],[4,5,-1],[1,2,4]]);
d1, u1 = eigen(rho1); sq_d1 = sqrt.(abs.(d1));
d2, u2 = eigen(rho2); sq_d2 = sqrt.(abs.(d2));
d3, u3 = eigen(rho3); sq_d3 = sqrt.(abs.(d3));
X1 = u1*diagm(0 => sq_d1)*u1'; X1inv = u1*diagm(0 => (1 ./sq_d1))*u1';
X2 = u2*diagm(0 => sq_d2)*u2'; X2inv = u2*diagm(0 => (1 ./sq_d2))*u2';
X3 = u3*diagm(0 => sq_d3)*u3'; X3inv = u3*diagm(0 => (1 ./sq_d3))*u3';
# execute gauge changes
Aprime = ncon(Any[A,X1,X2,X3],Any[[1,-2,2,3],[-1,1],[-3,2],[-4,3]]);
Bprime = ncon(Any[B,X1inv],Any[[-1,-2,1],[1,-3]]);
Fprime = ncon(Any[F,X2inv],Any[[-1,-2,1],[1,-3]]);
Cprime = ncon(Any[C,X3inv],Any[[-1,-2,1],[1,-3]]);
# new network is formed from tensors: {Aprime,Bprime,Cprime,D,E,Fprime,G}

# check both networks evaluate to the same tensor
connectlist = Any[[3,-5,4,5],[1,2,3],[6,-10,5],[-1,-2,1],[-3,-4,2],[-6,-7,4],[-8,-9,6]];
H0 = ncon(Any[A,B,C,D,E,F,G],connectlist)[:];
H1 = ncon(Any[Aprime,Bprime,Cprime,D,E,Fprime,G],connectlist)[:];
dH = norm(H0-H1) / norm(H0)