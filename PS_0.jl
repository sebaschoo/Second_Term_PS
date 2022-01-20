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

    return Π, y
    
end

vec_matrix = [zeros(3,3), 1:1:length(1), zeros(3,3), 1:1:length(1), zeros(3,3), 1:1:length(1), zeros(3,3), 1:1:length(1)]

for i in [0.2, 0.7, 0.9, 0.98]
    z = tauchen(3,i,0.4)
    if i==0.2
        vec_matrix[1] = z[1]
        vec_matrix[2] = z[2]
    elseif i==0.7
        vec_matrix[3] = z[1]
        vec_matrix[4] = z[2]
    elseif i==0.9
        vec_matrix[5] = z[1]
        vec_matrix[6] = z[2]
    else
        vec_matrix[7] = z[1]
        vec_matrix[8] = z[2]
    end
end  

y_tau = [zeros(1000,1), zeros(1000,1), zeros(1000,1), zeros(1000,1)]
pos = [zeros(1000,1), zeros(1000,1), zeros(1000,1), zeros(1000,1)]

for k=1:4

    l = 2*k - 1
    m = 2*k
    y_tau[k][1] = vec_matrix[m][rand(1:3)]

    for i=2:1000
        x=rand()
        for j=1:3
            if y_tau[k][i-1]==vec_matrix[m][j]
                pos[k][i] = j
            end
        end
        prob_row = vec_matrix[l][floor(Int, pos[k][i]),:]
        thr_1 = prob_row[1]
        thr_2 = thr_1 + prob_row[2]
        if x<=thr_1
            y_tau[k][i]=vec_matrix[m][1]
        elseif x>thr_1 && x<=thr_2
            y_tau[k][i] = vec_matrix[m][2]
        else
            y_tau[k][i] = vec_matrix[m][3]
        end
    end
end

t = collect(range(1,stop=1000,length=1000))
plot(t, y_tau[1])
plot(t, y_tau[2])
plot(t, y_tau[3])
plot(t, y_tau[4])

##############
# Question 2 #
##############

#Getting initial vector and probability matrix
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
ybar_rou = @show Vector(z[2])

#Simulating the time series
y_rou = zeros(1000,1)
pos = zeros(1000,1)
y_rou[1] = ybar_rou[rand(1:3)]

for i=2:1000
    x=rand()
    for j=1:3
        if y_rou[i-1]==ybar_rou[j]
            pos[i] = j
        end
    end
    prob_row = prob_rou[floor(Int, pos[i]),:]
    thr_1 = prob_row[1]
    thr_2 = prob_row[1] + prob_row[2]
    if x<=thr_1
        y_rou[i]=ybar_rou[1]
    elseif x>thr_1 && x<=thr_2
        y_rou[i] = ybar_rou[2]
    else
        y_rou[i] = ybar_rou[3]
    end
end

t = collect(range(1,stop=1000,length=1000))
plot(t, y_rou)

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
y_ar = normal_ar(ρ,N)
t = collect(range(1,stop=1000,length=1000))

plot(t, y_ar)


##############
# Question 4 #
##############

# Tauchen
z = tauchen(10,0.2,0.4)
prob_tau = z[1]
ybar_tau = @show Vector(z[2])

y_tau = zeros(1000,1)
pos = zeros(1000,1)
y_tau[1] = ybar_tau[rand(1:10)]

for i=2:1000
    x=rand()
    for j=1:10
        if y_tau[i-1]==ybar_tau[j]
            pos[i] = j
        end
    end
    prob_row = prob_tau[floor(Int, pos[i]),:]
    thr_1 = prob_row[1]
    thr_2 = thr_1 + prob_row[2]
    thr_3 = thr_2 + prob_row[3]
    thr_4 = thr_3 + prob_row[4]
    thr_5 = thr_4 + prob_row[5]
    thr_6 = thr_5 + prob_row[6]
    thr_7 = thr_6 + prob_row[7]
    thr_8 = thr_7 + prob_row[8]
    thr_9 = thr_8 + prob_row[9]
    if x<=thr_1
        y_tau[i]=ybar_tau[1]
    elseif x>thr_1 && x<=thr_2
        y_tau[i] = ybar_tau[2]
    elseif x>thr_2 && x<=thr_3
        y_tau[i] = ybar_tau[3]
    elseif x>thr_3 && x<=thr_4
        y_tau[i] = ybar_tau[4]
    elseif x>thr_4 && x<=thr_5
        y_tau[i] = ybar_tau[5]
    elseif x>thr_5 && x<=thr_6
        y_tau[i] = ybar_tau[6]
    elseif x>thr_6 && x<=thr_7
        y_tau[i] = ybar_tau[7]
    elseif x>thr_7 && x<=thr_8
        y_tau[i] = ybar_tau[8]
    elseif x>thr_8 && x<=thr_9
        y_tau[i] = ybar_tau[9]
    else
        y_tau[i] = ybar_tau[10]
    end
end

t = collect(range(1,stop=1000,length=1000))
plot(t, y_tau)

# Rouwenhorst
z = rouwenhorst(10,0.2,0.4,0)
prob_rou = z[1]
ybar_rou = @show Vector(z[2])

y_rou = zeros(1000,1)
pos = zeros(1000,1)
y_rou[1] = ybar_tau[rand(1:10)]

for i=2:1000
    x=rand()
    for j=1:10
        if y_rou[i-1]==ybar_rou[j]
            pos[i] = j
        end
    end
    prob_row = prob_rou[floor(Int, pos[i]),:]
    thr_1 = prob_row[1]
    thr_2 = thr_1 + prob_row[2]
    thr_3 = thr_2 + prob_row[3]
    thr_4 = thr_3 + prob_row[4]
    thr_5 = thr_4 + prob_row[5]
    thr_6 = thr_5 + prob_row[6]
    thr_7 = thr_6 + prob_row[7]
    thr_8 = thr_7 + prob_row[8]
    thr_9 = thr_8 + prob_row[9]
    if x<=thr_1
        y_rou[i]=ybar_rou[1]
    elseif x>thr_1 && x<=thr_2
        y_rou[i] = ybar_rou[2]
    elseif x>thr_2 && x<=thr_3
        y_rou[i] = ybar_rou[3]
    elseif x>thr_3 && x<=thr_4
        y_rou[i] = ybar_rou[4]
    elseif x>thr_4 && x<=thr_5
        y_rou[i] = ybar_rou[5]
    elseif x>thr_5 && x<=thr_6
        y_rou[i] = ybar_rou[6]
    elseif x>thr_6 && x<=thr_7
        y_rou[i] = ybar_rou[7]
    elseif x>thr_7 && x<=thr_8
        y_rou[i] = ybar_rou[8]
    elseif x>thr_8 && x<=thr_9
        y_rou[i] = ybar_rou[9]
    else
        y_rou[i] = ybar_rou[10]
    end
end

t = collect(range(1,stop=1000,length=1000))
plot(t, y_rou)


##############
# Question 5 #
##############

#= The coding for this question is already included previously in a loop
    =#
