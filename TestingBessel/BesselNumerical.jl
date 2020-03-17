import Pkg; using Pkg
Pkg.add("DifferentialEquations"); Pkg.add("Plots"); Pkg.add("SpecialFunctions"); Pkg.add("Calculus")
using DifferentialEquations, Plots, SpecialFunctions, Calculus


function besseleq(dv,v,y,p,t)
    a = p
    dv[1] = -v[1]/t - ( 1 - a^2/t^2 ) * y[1]  
    #Diffeq *** y'' = -y'/t - ( t^2 -a^2 ) * y
end

a=0.0
tspan=(1.0,10.0)
u0 = 1.0
du0 = 0.0
function solArray(somesol)
    z=zeros(length(somesol.t),1)
    length(z)
    for i in 1:1:length(z)
        z[i]=somesol.u[i][2]
    end
    return z 
end 

beq=SecondOrderODEProblem(besseleq,[du0], [u0], tspan, a)
sol=solve(beq,DPRKN12(),dense=true,abstol=1e-20,reltol=1e-20, dt=1.0, maxiters=1e50)
#plot(sol.t,solArray(sol))

nu=0.0
f(x)=besselj(nu,x)
function besselArray(somesol)
    z=zeros(length(somesol.t),1)
    length(z)
    for i in 1:1:length(z)
        z[i]=f(somesol.t[i]-1)-somesol.u[i][2] #exact - numerical
    end
    return z 
end
f(1.0009243471260376-1)
sol.u[2][2]
sol.t


using QuadGK
q=1.0027549075130426
BI(t)=cos(0*t-q*sin(t))/pi
quadgk(BI, 0, pi, rtol=1e-10)
f(q)
besselj(0,q)

#differentiate()