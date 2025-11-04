module MultifractalTools

using LinearAlgebra
using LsqFit
using Statistics



export obtain_qs, compute_scaling_quantities, compute_spectrum, plot_spectrum, plot_to_fit


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

function find_best_scaling_size(max_size::Integer, crop_ratio::Float64)

    min_search_size = floor(Int, max_size * (1.0 - crop_ratio))


    best_size = max_size 
    best_divs_list = get_divisors(best_size)
    best_count = length(best_divs_list)


    for n in min_search_size:max_size
        DivTemp = get_divisors(n)
        if length(DivTemp) >= best_count
            best_count = length(DivTemp)
            best_size = n
            best_divs_list = DivTemp
        end
    end

    return (size = best_size, divisors = best_divs_list)
end


function bin_data(data::AbstractMatrix{T}, l::Integer) where {T<:Number}

    
    SizeBinned = floor.(Int64,size(data) ./ l)
    BinnedData = zeros(T, SizeBinned...)

    for idx in CartesianIndices(data) 

        BoxIndex = (Tuple(idx) .- 1) .÷ l .+ 1
        BinnedData[BoxIndex...] += abs2(data[idx]) 
    end
    return BinnedData
end


function plot_to_fit()
    error("Plotting functionality requires GLMakie. Please run `using GLMakie` to load the plotting extension.")
end

function plot_spectrum()
    error("Plotting functionality requires GLMakie. Please run `using GLMakie` to load the plotting extension.")
end

function compute_scaling_quantities(
    data::AbstractMatrix{T},
    qs::AbstractVector{<:Real};
    ls::AbstractVector{<:Integer} = Integer[],
    crop_to_best_fit::Bool = true,
    crop_ratio::Float64 = 0.1
    ) where {T<:Number} 
    #The advanced user can provide specific ls values, otherwise we always compute them!



    

    #Define the sizes!
    best_size = minimum(size(data)) #Defining first on the outside scope 
    if isempty(ls)

        if crop_to_best_fit


        #Obtain the possible box sizes
        biggest_size_info = find_best_scaling_size(minimum(size(data)), crop_ratio)
        ls = biggest_size_info.divisors
        best_size = biggest_size_info.size
        else
            ls = get_divisors(minimum(size(data)))
        end
    
    end

    #Crop the data! 
    data_cropped = data[1:best_size, 1:best_size] #not centered, should be ok!


    Base.@info "Data cropped to $best_size x $best_size (found $(length(ls)) divisors)."
    #1. Renormalize the data so that sum of |data|^2 = 1.
    data_renorm = renormalize_data(data_cropped) 
    
    Zqs = zeros(T, length(ls), length(qs))
    S_qs = zeros(T, length(ls), length(qs))
    ZPrime_qs = zeros(T, length(ls), length(qs))





    for i_l in eachindex(ls)
        l = ls[i_l]
        binnedData = bin_data(data_renorm, l)
        miu = similar(binnedData)

        #Possibly add a filter to remogve 0 values - Later
        for i_q in eachindex(qs)
            q = qs[i_q]
            Zqs[i_l, i_q] = sum(binnedData .^ q)


            miu .= (binnedData .^ q) ./ Zqs[i_l, i_q]
            S_qs[i_l, i_q] = sum(miu .* log.(miu))
            ZPrime_qs[i_l,i_q] = sum(miu .* log.(binnedData))


        end


    end

    return (ls= ls , Zqs = Zqs, Sqs = S_qs, ZPrimes = ZPrime_qs ) #NamedTuple is complete amazing here!! Could be called as Tuple.Zqs etc
end


function power_law_model(x, p) #Log Scale
    return p[1] .* x .+ p[2]
end

function compute_spectrum(ScalingQuantities::NamedTuple, qs::AbstractVector{<:Real}, λ1::Integer, λ2::Integer)

    τqs = zeros(eltype(qs), length(qs))
    αqs = zeros(eltype(qs), length(qs))
    fqs = zeros(eltype(qs), length(qs))

    λs = log.(ScalingQuantities.ls ./ maximum(ScalingQuantities.ls))

    for i_q in eachindex(qs)


        logZs = log.(ScalingQuantities.Zqs[:, i_q])
        Ss = ScalingQuantities.Sqs[:, i_q]
        ZPrimes = ScalingQuantities.ZPrimes[:, i_q]  


        #Fit τ(q)
        p0 = [1.0,1.0]
        fit_τ = curve_fit(power_law_model, λs[λ1:λ2], logZs[λ1:λ2], p0) #I need to work on how to define this lambda1 or lambda2? Indices, Values?
        τqs[i_q] = fit_τ.param[1]

        #Fit α(q)  
        fit_α = curve_fit(power_law_model, λs[λ1:λ2], ZPrimes[λ1:λ2], p0)
        αqs[i_q] = fit_α.param[1]

        #Fit f(α)
        fit_f = curve_fit(power_law_model, λs[λ1:λ2], Ss[λ1:λ2], p0)
        fqs[i_q] = fit_f.param[1]

    end
    return (qs=qs, τqs = τqs, αs = αqs, fs = fqs)
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

end #module