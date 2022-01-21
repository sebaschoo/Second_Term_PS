                                            ################################
                                            ### Problem Set 0 - ECON 603 ###
                                            ################################

                                                #Sebastian Melo-Martin


# Call some useful packages 
using Distributions
using Random
using QuantEcon
using Plots
using NLopt
using SpecialFunctions
using StatsBase
using DataFrames
using DelimitedFiles
using CSV

Random.seed!(641993)

############
# Question 1
############

std_norm_cdf(x::T) where {T <: Real} = 0.5 * erfc(-x/sqrt(2))
std_norm_cdf(x::Array{T}) where {T <: Real} = 0.5 .* erfc(-x./sqrt(2))

function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::T3=3) where {T1 <: Real, T2 <: Real, T3 <: Real}
    # Discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # Fill in columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    yy = y .+ μ / (1 - ρ)  
    Π = Π./sum(Π, dims = 2)

    return Π, yy
    
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
plot(t, y_tau[1], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau3_1.png")
plot(t, y_tau[2], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau3_2.png")
plot(t, y_tau[3], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau3_3.png")
plot(t, y_tau[4], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau3_4.png")

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

vec_matrix = [zeros(3,3), 1:1:length(1), zeros(3,3), 1:1:length(1), zeros(3,3), 1:1:length(1), zeros(3,3), 1:1:length(1)]

for i in [0.2, 0.7, 0.9, 0.98]
    z = rouwenhorst(3,i,0.4,0)
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

y_rou = [zeros(1000,1), zeros(1000,1), zeros(1000,1), zeros(1000,1)]
pos = [zeros(1000,1), zeros(1000,1), zeros(1000,1), zeros(1000,1)]

#Simulating the time series
for k=1:4

    l = 2*k - 1
    m = 2*k
    y_rou[k][1] = vec_matrix[m][rand(1:3)]

    for i=2:1000
        x=rand()
        for j=1:3
            if y_rou[k][i-1]==vec_matrix[m][j]
                pos[k][i] = j
            end
        end
        prob_row = vec_matrix[l][floor(Int, pos[k][i]),:]
        thr_1 = prob_row[1]
        thr_2 = thr_1 + prob_row[2]
        if x<=thr_1
            y_rou[k][i]=vec_matrix[m][1]
        elseif x>thr_1 && x<=thr_2
            y_rou[k][i] = vec_matrix[m][2]
        else
            y_rou[k][i] = vec_matrix[m][3]
        end
    end
end

t = collect(range(1,stop=1000,length=1000))
plot(t, y_rou[1], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou3_1.png")
plot(t, y_rou[2], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou3_2.png")
plot(t, y_rou[3], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou3_3.png")
plot(t, y_rou[4], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou3_4.png")


##############
# Question 3 #
##############

function normal_ar(ρ, N)
    d = Normal(0, 0.16)
    y = [0.0 for i = 1:N]
    ϵ = rand(d, N)
    for i in 1:(N-1)
        y[i+1] = ρ*y[i] + ϵ[i] 
    end
    return y
end

y_ar = [zeros(1000), zeros(1000), zeros(1000), zeros(1000)]

N=1000
for i in [0.2, 0.7, 0.9, 0.98]
    ρ=i
    if i==0.2
        y_ar[1] = normal_ar(ρ,N)
    elseif i==0.7
        y_ar[2] = normal_ar(ρ,N)
    elseif i==0.9
        y_ar[3] = normal_ar(ρ,N)
    else
        y_ar[4] = normal_ar(ρ,N)
    end 
end

t = collect(range(1,stop=1000,length=1000))
plot(t, y_ar[1], title=" ", legend = :false, xlabel = "Time", color = :tan2)
savefig("./PS0_ar_1.png")
plot(t, y_ar[2], title=" ", legend = :false, xlabel = "Time", color = :tan2)
savefig("./PS0_ar_2.png")
plot(t, y_ar[3], title=" ", legend = :false, xlabel = "Time", color = :tan2)
savefig("./PS0_ar_3.png")
plot(t, y_ar[4], title=" ", legend = :false, xlabel = "Time", color = :tan2)
savefig("./PS0_ar_4.png")


##############
# Question 4 #
##############

# Tauchen
vec_matrix = [zeros(10,10), zeros(10), zeros(10,10), zeros(10), zeros(10,10), zeros(10), zeros(10,10), zeros(10)]

for i in [0.2, 0.7, 0.9, 0.98]
    z = tauchen(10,i,0.4)
    if i==0.2
        vec_matrix[1] = z[1]
        vec_matrix[2] = Vector(z[2])
    elseif i==0.7
        vec_matrix[3] = z[1]
        vec_matrix[4] = Vector(z[2])
    elseif i==0.9
        vec_matrix[5] = z[1]
        vec_matrix[6] = Vector(z[2])
    else
        vec_matrix[7] = z[1]
        vec_matrix[8] = Vector(z[2])
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
        for j=1:10
            if y_tau[k][i-1]==vec_matrix[m][j]
                pos[k][i] = j
            end
        end
        prob_row = vec_matrix[l][floor(Int, pos[k][i]),:]
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
            y_tau[k][i]= vec_matrix[m][1]
        elseif x>thr_1 && x<=thr_2
            y_tau[k][i] = vec_matrix[m][2]
        elseif x>thr_2 && x<=thr_3
            y_tau[k][i] = vec_matrix[m][3]
        elseif x>thr_3 && x<=thr_4
            y_tau[k][i] = vec_matrix[m][4]
        elseif x>thr_4 && x<=thr_5
            y_tau[k][i] = vec_matrix[m][5]
        elseif x>thr_5 && x<=thr_6
            y_tau[k][i] = vec_matrix[m][6]
        elseif x>thr_6 && x<=thr_7
            y_tau[k][i] = vec_matrix[m][7]
        elseif x>thr_7 && x<=thr_8
            y_tau[k][i] = vec_matrix[m][8]
        elseif x>thr_8 && x<=thr_9
            y_tau[k][i] = vec_matrix[m][9]
        else
            y_tau[k][i] = vec_matrix[m][10]
        end
    end
end   

t = collect(range(1,stop=1000,length=1000))
plot(t, y_tau[1], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau10_1.png")
plot(t, y_tau[2], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau10_2.png")
plot(t, y_tau[3], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau10_3.png")
plot(t, y_tau[4], title=" ", legend = :false, xlabel = "Time")
savefig("./PS0_tau10_4.png")

# Rouwenhorst
vec_matrix = [zeros(10,10), zeros(10), zeros(10,10), zeros(10), zeros(10,10), zeros(10), zeros(10,10), zeros(10)]

for i in [0.2, 0.7, 0.9, 0.98]
    z = rouwenhorst(10,i,0.4,0)
    if i==0.2
        vec_matrix[1] = z[1]
        vec_matrix[2] = Vector(z[2])
    elseif i==0.7
        vec_matrix[3] = z[1]
        vec_matrix[4] = Vector(z[2])
    elseif i==0.9
        vec_matrix[5] = z[1]
        vec_matrix[6] = Vector(z[2])
    else
        vec_matrix[7] = z[1]
        vec_matrix[8] = Vector(z[2])
    end
end  

y_rou = [zeros(1000,1), zeros(1000,1), zeros(1000,1), zeros(1000,1)]
pos = [zeros(1000,1), zeros(1000,1), zeros(1000,1), zeros(1000,1)]


for k=1:4

    l = 2*k - 1
    m = 2*k
    y_rou[k][1] = vec_matrix[m][rand(1:3)]

    for i=2:1000
        x=rand()
        for j=1:10
            if y_rou[k][i-1]==vec_matrix[m][j]
                pos[k][i] = j
            end
        end
        prob_row = vec_matrix[l][floor(Int, pos[k][i]),:]
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
            y_rou[k][i]=vec_matrix[m][1]
        elseif x>thr_1 && x<=thr_2
            y_rou[k][i] = vec_matrix[m][2]
        elseif x>thr_2 && x<=thr_3
            y_rou[k][i] = vec_matrix[m][3]
        elseif x>thr_3 && x<=thr_4
            y_rou[k][i] = vec_matrix[m][4]
        elseif x>thr_4 && x<=thr_5
            y_rou[k][i] = vec_matrix[m][5]
        elseif x>thr_5 && x<=thr_6
            y_rou[k][i] = vec_matrix[m][6]
        elseif x>thr_6 && x<=thr_7
            y_rou[k][i] = vec_matrix[m][7]
        elseif x>thr_7 && x<=thr_8
            y_rou[k][i] = vec_matrix[m][8]
        elseif x>thr_8 && x<=thr_9
            y_rou[k][i] = vec_matrix[m][9]
        else
            y_rou[k][i] = vec_matrix[m][10]
        end
    end
end

t = collect(range(1,stop=1000,length=1000))
plot(t, y_rou[1], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou10_1.png")
plot(t, y_rou[2], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou10_2.png")
plot(t, y_rou[3], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou10_3.png")
plot(t, y_rou[4], title=" ", legend = :false, xlabel = "Time", color = :deepskyblue3)
savefig("./PS0_rou10_4.png")

table1 = zeros(4,6)
table2 = zeros(4,6)
table3 = zeros(4,6)

for j=1:4
    table1[j,1] = mean(y_tau[j])
    table1[j,2] = var(y_tau[j])
    table1[j,3] = percentile(vec(y_tau[j]),25)
    table1[j,4] = percentile(vec(y_tau[j]),50)
    table1[j,5] = percentile(vec(y_tau[j]),75)
    table1[j,6] = percentile(vec(y_tau[j]),90)
end 

CSV.write("table1.csv",  Tables.table(table1), writeheader=false)

for j=1:4
    table2[j,1] = mean(y_rou[j])
    table2[j,2] = var(y_rou[j])
    table2[j,3] = percentile(vec(y_rou[j]),25)
    table2[j,4] = percentile(vec(y_rou[j]),50)
    table2[j,5] = percentile(vec(y_rou[j]),75)
    table2[j,6] = percentile(vec(y_rou[j]),90)
end 

CSV.write("table2.csv",  Tables.table(table2), writeheader=false)

for j=1:4
    table3[j,1] = mean(y_ar[j])
    table3[j,2] = var(y_ar[j])
    table3[j,3] = percentile(vec(y_ar[j]),25)
    table3[j,4] = percentile(vec(y_ar[j]),50)
    table3[j,5] = percentile(vec(y_ar[j]),75)
    table3[j,6] = percentile(vec(y_ar[j]),90)
end 

CSV.write("table3.csv",  Tables.table(table3), writeheader=false)

##############
# Question 5 #
##############

#= The coding for this question is already included previously in loops
    =#