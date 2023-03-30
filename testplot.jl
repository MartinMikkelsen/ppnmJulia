using FewSpecialFunctions

x = range(0,500,500)

plot(x,FermiDiracIntegralNorm.(3/2,x))