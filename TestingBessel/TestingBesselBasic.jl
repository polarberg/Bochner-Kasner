import Pkg; using Pkg
Pkg.add("SpecialFunctions"); Pkg.add("Plots"); #Pkg.add("MTH229");
using SpecialFunctions, Plots

#Bessel Fx of order 1
nu=.5
f(x)=besselj(nu,x)
g(x)=besselj(nu+1,x)
k(x)=besselj(nu+2,x)
plot([f,g,k],-10,10)

#Exact Sol vs Calculated Sol of J_0.5
J_half(x)=sqrt(2/(pi*x))*sin(x)
y=1
e(y)=besselj(0.5,y)
n(y)=J_half(y)
d(y)=e(y)-n(y)
plot(d,0,10)



sphericalbesselj(ν, x) = √(π/2x)*besselj(ν+1/2, x)
sphericalbessely(ν, x) = √(π/2x)*bessely(ν+1/2, x)



start=  0
finish =   10
span=[]
plot()
for n in Int32[0,2]
    print(typeof(n))
    b(x)=besselj(n,x)
    plot(b,0,10)
    plot(b,start,finish)
end