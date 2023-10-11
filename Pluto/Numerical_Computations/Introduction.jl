
@show p = 22/7
@show float(π)
acc = abs(p-π)
println("absolute accuracy = $acc")
println("relative accuracy = $(acc/π)")
println("Number of accurate digits = $(-log10(acc/π))")
@show eps()
@show eps(161.8)

"""
    horner(c,x)

Evaluate a polynomial whose coefficients are given in ascending
order in `c`, at the point `x`, using Horner's rule.
"""
function horner(c,x)
    n = length(c)
    y = c[n]
    for k in n-1:-1:1
        y = x*y + c[k]
    end
    return y
end

year = [1982,2000,2010,2015]; 
pop = [1008.18, 1262.64, 1337.82, 1374.62];
t = year .- 1980.0
y = pop;
V = [ t[i]^j for i=1:4, j=0:3 ]
c = V \ y
y - V*c