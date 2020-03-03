#DEPOT_PATH
#DEPOT_PATH[1]="D:\\.julia"
#DEPOT_PATH[2]="D:\\Julia-1.2.0\\local\\share\\julia"
#DEPOT_PATH[3]="D:\\Julia-1.2.0\\share\\julia"
#DEPOT_PATH

import Pkg; using Pkg
Pkg.add("DifferentialEquations"); Pkg.add("Plots")
using DifferentialEquations, Plots #http://docs.juliadiffeq.org/latest/ #http://docs.juliaplots.org/latest/

function waveeq(dv,v,y,p,t)
    kx,kz = p
    dv[1] = v[1]/(3*t) - ( (kx^2)/t^(4/3) + kz^2*t^(2/3) ) * y[1]  #Diffeq *** y'' = y'/(3t) - ( (kx^2)/t^(4/3) + kz^2 * t^2*t^(2/3) ) * y
    #dv[1] = -y[1]
end

kx=.7;kz=1.1; kx = big.(kx); kz=big.(kz)
    tspan = (4.0,1e-16); dt=-1
    #tspan = (4.0,12); 
    u0range=range(-1,1,step=.2)
    #a=Vector(u0range)
    #b=Vector(u0range)
    a=[0.0,1.0,0.5,0.5]
    b=[1.0,0.0,0.5,-0.5]
u0 = big.(a)
du0 = big.(b)
function solArray(somesol)
    z=zeros(length(somesol.t),1)
    length(z)
    for i in 1:1:length(z)
        z[i]=somesol.u[i][2]
    end
    return z 
end 
#
function solveODE(n,i)    
    solve(SecondOrderODEProblem(waveeq, [du0[n]], [u0[i]], big.(tspan), (big.(kx),big.(kz))), DPRKN6(), dt=-1,abstol=1e-10,reltol=1e-10, tstops=0.0) # 
end
sol = solveODE(2,2)


sol1=solve(SecondOrderODEProblem(waveeq, [1.0], [0.0], big.(tspan), (big.(kx),big.(kz))), DPRKN6(),dt=-1, abstol=1e-10,reltol=1e-10, tstops=0.0) # 
plot!(sol1.t,solArray(sol1),title="testing one case")
sol2=solve(SecondOrderODEProblem(waveeq, [0.0], [0.65], big.(tspan), (big.(kx),big.(kz))), DPRKN6(),dt=-1, abstol=1e-10,reltol=1e-10, tstops=0.0) # 
plot!(sol2.t,solArray(sol2),title="testing one case")
sol3=solve(SecondOrderODEProblem(waveeq, [0.5], [0.5], big.(tspan), (big.(kx),big.(kz))), DPRKN6(),dt=-1, abstol=1e-10,reltol=1e-10, tstops=0.0) # 
plot!(sol3.t,solArray(sol3),title="testing one case")
sol4=solve(SecondOrderODEProblem(waveeq, [-0.5], [0.5], big.(tspan), (big.(kx),big.(kz))), DPRKN6(),dt=-1, abstol=1e-10,reltol=1e-10, tstops=0.0) # 
plot!(sol4.t,solArray(sol4),title="testing one case")



#solArray(sol)



function PlotAppend(p::Plots.Plot,NumofPlots=length(a))	
    for n in 1:length(a)
        for i in 1:length(a)
            tdu0=fill(a[n],length(a))           
            tsol=solveODE(n,i)
            plot!(tsol.t,solArray(tsol),legend = false)
            #plot!(solveODE(n,i),legend = false) # , label=string("Numerical sol ",n)    ,ls=:dashdot
        end
    end
    
    return plots
end
plots = plot(title="Kasner 2/3, 2/3, -1/3 Solution, t=1e-16 to t=10, u0 at t=2", xlabel="time")  #Construct base object
PlotAppend(plots)

# saving the plot 
plot!(sol1.t,solArray(sol1),legend = false,dpi=2000)
#png("D:\\COMP SCI\\Julia\\Final\\Finalfigure2.png")  #**** CHANGE to WHERE YOU want to SAVE
#savefig("C:\\Users\\Austin\\Documents\\CS\\Julia\\Kasner\\Kasner_2-3_2-3_-1-3__1e-16to12__u0_at_t-4.pdf") #*** remember the // for each / 