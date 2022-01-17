################################
### Problem Set 0 - ECON 603 ###
################################

#Name: Sebastian Melo-Martin

#1

#2

#3

function normal_ar(ρ, N)

    d = Normal(0, 0.16)

    y = [0.0 for i = 1:N]

    ϵ = rand(d, N)

    for i in 1:(N-1)
        y[i+1] = ρ*y[i] + ϵ[i] 
    end

    return y
end

ρ=0.2
N=1000
z = normal_ar(ρ,N)

plot!(1:1000, z)
z