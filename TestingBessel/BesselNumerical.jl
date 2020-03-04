import Pkg; using Pkg
Pkg.add("DifferentialEquations"); Pkg.add("Plots"); Pkg.add("Calculus")
using DifferentialEquations, Plots, SpecialFunctions, Calculus


function besseleq(dv,v,y,p,t)
    a = p
    dv[1] = -v[1]/t - ( 1 - a^2/t^2 ) * y[1]  
    #Diffeq *** y'' = -y'/t - ( t^2 -a^2 ) * y
end

nu=1
besselj(1,1)
bessely(1,1)

differentiate()



a=1.0
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

sol=solve(SecondOrderODEProblem(besseleq,[du0],[u0],tspan,(a)),DP)
plot(sol.t,solArray(sol))
