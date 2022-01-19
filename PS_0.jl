################################
### Problem Set 0 - ECON 603 ###
################################

#Name: Sebastian Melo-Martin

using Distributions
using Random
using QuantEcon

#1

function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::T3=3) where {T1 <: Real, T2 <: Real, T3 <: Real}
    # Discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = cdf(Normal(),((y[1] - ρ*y[row] + d/2) / σ))
        Π[row, N] = 1 - cdf(Normal(),((y[N] - ρ*y[row] - d/2) / σ))

        # Fill in columns
        for col = 2:N-1
            Π[row, col] = (cdf(Normal(),((y[col] - ρ*y[row] + d/2) / σ)) -
                           cdf(Normal(),((y[col] - ρ*y[row] - d/2) / σ)))
        end
    end

    MarkovChain(Π, y)
end

tauchen(2,0.2,0.4)
z

#2

#3

function normal_ar(ρ, N)
    d = Normal(0, 0.16)
    y = [0.0 for i = 1:N]
    Random.seed!(641993)
    ϵ = rand(d, N)
    for i in 1:(N-1)
        y[i+1] = ρ*y[i] + ϵ[i] 
    end
    return y
end

ρ=0.2
N=1000
z = normal_ar(ρ,N)
t = collect(range(1,stop=1000,length=1000))

plot(t, z)
