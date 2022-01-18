################################
### Problem Set 0 - ECON 603 ###
################################

#Name: Sebastian Melo-Martin

#1

function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::T3=3) where {T1 <: Real, T2 <: Real, T3 <: Real}
    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # fill in the middle columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    yy = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable

    # renormalize. In some test cases the rows sum to something that is 2e-15
    # away from 1.0, which caused problems in the MarkovChain constructor
    Π = Π./sum(Π, dims = 2)

    MarkovChain(Π, yy)
end

z = tauchen(1000,0.2,0.4)


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
