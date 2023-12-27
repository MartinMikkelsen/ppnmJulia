using LinearAlgebra, TensorOperations

d = 20;
A = rand(d,d,d); B = rand(d,d,d); C = rand(d,d,d);

##### Evaluate network via index summation
function tempfunct(A,B,C,d)
    D0 = zeros(d,d,d);
    for b1 = 1:d
        for a2 = 1:d
            for c3 = 1:d
                for a1 = 1:d
                    for a3 = 1:d
                        for c1 = 1:d
                            D0[b1,a2,c3] = D0[b1,a2,c3]+A[a1,a2,a3]*B[b1,a1,c1]*C[c1,a3,c3];
                        end
                    end
                end
            end
        end
    end
    return D0
end

t_sum = @elapsed D0 = tempfunct(A,B,C,d);

##### Evaluate network using reshape and permute
function tempfunct2(A,B,C,d)
    Xmid = reshape(reshape(permutedims(B,[1,3,2]),d^2,d)*reshape(A,d,d^2),d,d,d,d);
    D1 = reshape(reshape(permutedims(Xmid,[1,3,2,4]),d^2,d^2)*reshape(C,d^2,d),d,d,d);
    return D1
end
t_res = @elapsed D1 = tempfunct2(A,B,C,d)

##### Evaluate using ncon
t_ncon = @elapsed D2 = ncon(Any[A,B,C],Any[[1,-2,2],[-1,1,3],[3,2,-3]]; order = [1,2,3]);

##### Compare
tdiffs = [maximum(abs.(D0[:]-D1[:])),maximum(abs.(D1[:]-D2[:])),maximum(abs.(D2[:]-D0[:]))]
ttimes = [t_sum, t_res, t_ncon]