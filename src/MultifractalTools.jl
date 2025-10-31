module MultifractalTools

using LinearAlgebra
using LsqFit
using Statistics


export obtain_qs, compute_scaling_quantities


function renormalize_data(data::AbstractMatrix{T}) where T<:Number
    value = sum(abs2.(data))
    return data / sqrt(value)
end

function get_divisors(n::Integer) 
    n = abs(n)
    divs = Int[]
    for i in 1:floor(Int, √n)
        if n % i == 0
            push!(divs, i)
            if i != n ÷ i
                push!(divs, n ÷ i)
            end
        end
    end
    sort!(divs)
    return divs
end


# Find the number up to `cutoff` with the most divisors
function get_biggest_divisor(cutoff::Integer) 
    best_number = 0
    best_count = 0
    Div_List = []
    for n in 1:cutoff
        DivTemp = get_divisors(n)
        if length(DivTemp) >= best_count
            best_count = length(DivTemp)
            best_number = n
            if !isempty(Div_List)
                deleteat!(Div_List, 1)
            end
            push!(Div_List, DivTemp)
        end
    end
    return best_number, Div_List[end]
end


function bin_data(data::AbstractMatrix{T}, l::Integer) where {T<:Number}

    data = RenormalizeData(data) #The user should not do this, we do it here to be sure!

    SizeBinned = floor.(Int64,size(data) ./ l)
    BinnedData = zeros(T, SizeBinned...)

    for idx in CartesianIndices(data) 

        BoxIndex = (Tuple(idx) .- 1) .÷ l .+ 1
        BinnedData[BoxIndex...] += abs2(data[idx]) 
    end
    return BinnedData
end


function compute_scaling_quantities(data::AbstractMatrix{T},qs::AbstractVector{<:Real}, ls::AbstractVector{<:Integer}) where {T<:Number} 

    #1. Renormalize the data so that sum of |data|^2 = 1.
    data_renorm = renormalize_data(data) 

    #2. Define the output matrix

    #ls should be computed here!! The user should not have to think on this (only if he wants - kwarg maybe?)
    Zqs = zeros(T, length(ls), length(qs))
    S_qs = zeros(T, length(ls), length(qs))
    Zprime_qs = zeros(T, length(ls), length(qs))





    for i_l in eachindex(ls)
        l = ls[i_l]
        binnedData = bin_data(data_renorm, l)
        miu = similar(T, size(binnedData))        

        #Possibly add a filter to remogve 0 values - Later
        for i_q in eachindex(qs)
            q = qs[i_q]
            Zqs[i_l, i_q] = sum(binnedData .^ q)


            miu .= (binnedData .^ q) ./ Zqs[i_l, i_q]
            S_qs[i_l, i_q] = sum(miu .* log.(miu))
            ZPrime_qs[i_l,i_q] = sum(miu .* log.(binnedData))


        end


    end

    return (Zqs = Zqs, Sqs = S_qs, ZPrimes = ZPrime_qs ) #NamedTuple is complete amazing here!! Could be called as Tuple.Zqs etc
    
    return 
end


function power_law_model(x, p) #Log Scale
    return p[1] .* x .+ p[2]
end

function compute_spectrum(ScalingQuantities::NamedTuple, qs::AbstractVector{<:Real}, ls::AbstractVector{<:Integer}, λ1::Real, λ2::Real)

    τqs = zeros(eltype(qs), length(qs))
    αqs = zeros(eltype(qs), length(qs))
    fqs = zeros(eltype(qs), length(qs))


    for i_q in eachindex(qs)
        q = qs[i_q]

        λs = log.(ls ./ maximum(ls))
        Zs = ScalingQuantities.Zqs[:, i_q]
        Ss = ScalingQuantities.Sqs[:, i_q]
        logZPrimes = ScalingQuantities.ZPrimes[:, i_q]  


        #Fit τ(q)
        p0 = [1.0,1.0]
        fit_τ = curve_fit(power_law_model, log.(λs[λ1:λ2]), logZs[λ1:λ2], p0) #I need to work on how to define this lambda1 or lambda2? Indices, Values?
        τqs[i_q] = fit_τ.param[1]

        #Fit α(q)  
        fit_α = curve_fit(power_law_model, log.(λs[λ1:λ2]), Ss[λ1:λ2], p0)
        αqs[i_q] = fit_α.param[1]

        #Fit f(α)
        fit_f = curve_fit(power_law_model, log.(λs[λ1:λ2]), logZPrimes[λ1:λ2], p0)
        fqs[i_q] = fit_f.param[1]

    end
    return (τqs = τqs, αs = αqs, fs = fqs)
end




function obtain_qs(qmin::Number, qmax::Number, num_q::Integer) 
    return collect(LinRange(qmin, qmax, num_q)) #This could be memory-problematic but for now it's ok
end



#= 
For the user, the workflow should be:

1) Choose the q values you want to work with 

qs = ObtainQs(qmin, qmax, num_q) #Done

2) Obtain the curves of the partition function for all box sizes. 

Zqs = ComputeAllParition(data, qs)


3) Choose the λ_1 and λ_2 values to perform the fits

I need to think better on how to do this automatically 
    - manual approach is always the safest, but how can we handle large datasets? 
    - Some minimization of the error? 
    - Can this lead to wrong results (catching local minima basicaly)?

4) Obtain tau(q) and f(α) 
τ, α, f = ObtainMultifractalSpectra(Zqs, qs, λ1, λ2)
=#

end