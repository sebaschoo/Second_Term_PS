                                            ################################
                                            ### Problem Set 0 - ECON 603 ###
                                            ################################

                                                #Sebastian Melo-Martin


# Call some useful packages 
using Distributions
using Random
using QuantEcon
using Plots


############
# Question 1
############

function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::T3=4) where {T1 <: Real, T2 <: Real, T3 <: Real}
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

    return Π, y
    
end

z = tauchen(3,0.2,0.4,3)
prob_tau = z[1]
y_tau = @show Vector(z[2])

y = zeros(1000,1)
pos = zeros(1000,1)
y[1] = y_tau[rand(1:3)]
y

for i=2:1000

    x=rand()
    for j=1:3
        if y[i-1]==y_tau[j]
            pos[i] = j
        end
    end
    prob_row = prob_tau[floor(Int, pos[i]),:]
    thr_1 = prob_row[1]
    thr_2 = thr_1 + prob_row[1]
    if x<=thr_1
        y[i]=y_tau[1]
    elseif x>thr_1 && x<=thr_2
        y[i] = y_tau[2]
    else
        y[i] = y_tau[3]
    end

end

t = collect(range(1,stop=1000,length=1000))
plot(t, y)

##############
# Question 2 #
##############

function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real)

    σ_y = σ / sqrt(1-ρ^2)
    p  = (1+ρ)/2
    Θ = [p 1-p; 1-p p]
    ψ = sqrt(N-1) * σ_y
    m = μ / (1 - ρ)

    state_values, p = _rouwenhorst(p, p, m, ψ, N)

    return p, state_values
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m-Δ, m+Δ],  [p 1-p; 1-q q]
    else
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n-1)
        θN = p    *[θ_nm1 zeros(n-1, 1); zeros(1, n)] +
             (1-p)*[zeros(n-1, 1) θ_nm1; zeros(1, n)] +
             q    *[zeros(1, n); zeros(n-1, 1) θ_nm1] +
             (1-q)*[zeros(1, n); θ_nm1 zeros(n-1, 1)]

        θN[2:end-1, :] ./= 2

        return range(m-Δ, stop=m+Δ, length=n), θN
    end
end

z = rouwenhorst(3,0.2,0.4,0)
prob_rou = z[1]

y_rou = @show Vector(z[2])

##############
# Question 3 #
##############

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


##############
# Question 4 #
##############

# Tauchen

# Rouwenhorst

# AR(1)
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
N=100
z = normal_ar(ρ,N)
t = collect(range(1,stop=1000,length=1000))

plot(t, z)

##############
# Question 5 #
##############