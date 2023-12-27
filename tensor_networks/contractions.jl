using LinearAlgebra

A = rand(4,4)

B = Matrix{Float64}(I, 5, 5)

C = ones(2,4,2,4)

D = zeros(3,5)

E = rand(2,3,4) + im*rand(2,3,4)

AP = permutedims(A,[1,2])

Q = rand(4,4,4)

BR = reshape(Q, 4,4^2)

d = 10;

A = rand(d,d,d,d); B = rand(d,d,d,d);

Ap = permutedims(A,[1,3,2,4])
Bp = permutedims(B,[1,4,2,3])
App = reshape(Ap,d^2,d^2); Bpp = reshape(Bp,d^2,d^2);
Cpp = App*Bpp
C = reshape(Cpp,d,d,d,d)

### Series of binary contractions


d = 10; A = rand(d,d);  B = rand(d,d); C = rand(d,d);

F0 = zeros(d,d)
for id = 1:d
    for jd = 1:d
        for kd = 1:d
            for ld = 1:d
                F0[id,jd] = F0[id,jd] + A[id,kd]*B[kd,ld]*C[ld,jd];
            end
        end
    end
end

F0

F1 = (A*B)*C

### Using ncon

using TensorOperations

d = 10;
A = rand(d,d,d); B = rand(d,d,d,d);
C = rand(d,d,d); D = rand(d,d);

TensorArray = Any[A,B,C,D];
IndexArray = Any[[1,-2,2],[-1,1,3,4],[5,3,2],[4,5]]

E = ncon(TensorArray,IndexArray, order = [5,3,4,1,2])

### Note that cont_order has been renamed to order

A = rand(d,d,d,d,d,d); 

B = ncon(Any[A],Any[[-1,-2,1,-3,-4,1]])

A = rand(d,d);
B = rand(d,d);

C = ncon(Any[A,B],Any[[-1,-2],[-3,-4]])