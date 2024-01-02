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