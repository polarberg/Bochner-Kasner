#DEPOT_PATH
#DEPOT_PATH[1]="D:\\.julia"
#DEPOT_PATH[2]="D:\\Julia-1.2.0\\local\\share\\julia"
#DEPOT_PATH[3]="D:\\Julia-1.2.0\\share\\julia"
#DEPOT_PATH

import Pkg; using Pkg
Pkg.add("TaylorSeries")
using TaylorSeries


t = Taylor1(Float64, 5)
y=1/(1-t)
#t=exp(t)
convert(Taylor1{Rational{BigInt}},y)
show_monomials(3)




# https://www.theoj.org/joss-papers/joss.01043/10.21105.joss.01043.pdf
# https://julialang.org/blog/2017/03/piday calculating Pi

#=
The example above can be slightly modified to compute, for example, the 100th Hermite
polynomial. In this case, the coefficients will be larger than 2
63 − 1, so the modular
behavior, under overflow of the standard Int64 type, will not suffice. Rather, the
polynomials should be generated with hermite_polynomials(BigInt, 100) to ensure
the use of arbitrary-length integers.
=#

displayBigO(false)
function hermite_polynomials(::Type{T}, nmax::Int) where {T <: Integer}
    x = Taylor1(T, nmax) # Taylor variable
    H = fill(x, nmax + 1) # vector of Taylor series to be overwritten
    H[1] = 1 # order 0
    H[2] = 2x # order 1
    for n in 2:nmax
        # recursion relation for order n:
        H[n+1] = 2x * H[n] - 2(n-1) * H[n-1]
    end
    return H
end
hermite_polynomials(n) = hermite_polynomials(Int, n);
H = hermite_polynomials(10);
function hermite_polynomial(n::Int)
    @assert 0 <= n <= length(H) "Not enough Hermite polynomials generated"
    return H[n+1]
end
hermite_polynomial(6)




#=
As a second example, we describe a numerical way of obtaining the Hermite polynomials from their 
generating function: the nth Hermite polynomial corresponds to the nth
derivative of the function exp(2t x − t
2).

This example shows that the calculations are performed numerically and not symbolically, using TaylorSeries.jl as a polynomial manipulator; this is manifested by the
fact that the last coefficient of HH(6) is not identical to an integer.
=#

G(x,t) = exp(2t * x - t^2) # generating function; G is typed as \scrG<TAB>
xn = set_variables("x", numvars=1, order=10)
x = xn[1]
t = Taylor1([zero(x), one(x)], 10) # Taylor1{TaylorN{Float64}}
gf = G(x, t) # Taylor1 expansion of G
HH(n::Int) = derivative(n, gf) # n-th derivative of `gf`
HH(6)

#=
Taylor method for integrating ordinary differential equations
    As a final example, we give a simple implementation of Picard iteration to integrate an
    ordinary differential equation, which is equivalent to the Taylor method.
    We consider the initial-value problem x˙ = x, with initial condition x(0) = 1. One step
    of the integration corresponds to constructing the Taylor series of the solution x(t) in
    powers of t:
=#    

∫_dt(u::Taylor1) = integrate(u) # the symbol R is obtained as \int<TAB>
function taylor_step(f, u0)
    u = copy(u0)
    unew = u0 + ∫_dt(f(u0))
    
    while unew != u
        u = unew
        unew = u0 + ∫_dt(f(u)) # Picard iteration
    end

    return u
end

f(x) = -x # Differential equation
order = 20 # maximum order of the Taylor expansion for the solution
u0 = Taylor1([1.0], order) # initial condition given as a Taylor expansion
solution = taylor_step(f, u0); # solution
exactsol(x) = exp(-x)
t=1
solution
solution(1)
solution(t) - exactsol(t) # compare solution with the exact value at t=1

#=
Thus this Taylor expansion of order 20 around t0 = 0 suffices to obtain the exact solution
at t = 1, while the error at time t = 2 from the same expansion is 4.53 × 10−14. This
indicates that a proper treatment should estimate the size of the required step that
should be taken as a function of the solution.
=#