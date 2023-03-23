using QuasiMonteCarlo
m = 4
d = 3
b = QuasiMonteCarlo.nextprime(d)
N = b^m # Number of points
pad = m # Length of the b-ary decomposition number = sum(y[k]/b^k for k in 1:pad)

# Unrandomized low discrepency sequence
x_faure = QuasiMonteCarlo.sample(N, d, FaureSample())

# Randomized version
x_nus = randomize(x_faure, QuasiMonteCarlo.OwenScramble(base = b, pad = pad)) # equivalent to sample(N, d, FaureSample(R = OwenScramble(base = b, pad = pad)))
x_lms = randomize(x_faure, MatousekScramble(base = b, pad = pad))
x_digital_shift = randomize(x_faure, DigitalShift(base = b, pad = pad))
x_shift = randomize(x_faure, Shift())
x_uniform = rand(d, N) # plain i.i.d. uniform