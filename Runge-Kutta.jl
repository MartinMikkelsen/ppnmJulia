using Plots, BenchmarkTools

function euler(dydt, y0, t; args=())
    ndof = length(y0)
    ntimes = length(t)
    y = zeros(ndof, ntimes)
    y[:, 1] = y0
    for cont = 2:ntimes
        h = t[cont] - t[cont - 1]
        y[:, cont] = y[:, cont - 1] + h*dydt(y[:, cont - 1], t[cont], args...)
    end
    return y
end

function pend(y, t, b, c)
    theta, omega = y
    dydt = [omega, -b*omega - c*sin(theta)]
    return dydt
end


b = 0.25
c = 5.0
y0 = [pi - 0.1, 0.0]
times = range(0, 10, 1001)
sol = euler(pend, y0, times, args=(b, c))

plot(times,sol[1,:],label="ω(t)",linewidth=2.5)
plot!(times,sol[2,:],label="θ(t)",linewidth=2.5)

@benchmark euler(pend, y0, times, args=(b, c))
