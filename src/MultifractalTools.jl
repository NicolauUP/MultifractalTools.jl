module MultifractalTools

using LinearAlgebra
using LsqFit
using Statistics


export nothing


function RenormalizeData(data::AbstractMatrix{T}) where T<:Number
    value = sum(abs2.(data))
    return data / sqrt(value)
end

function GetDivisors(n::Number{T}) where T<:Integer
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
function GetBiggestDivider(cutoff::Number{T}) where T<:Integer
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


function BinData(data::AbstractMatrix{T}, l::Number{V}) where {T<:Number,V<:Integer}

    data = RenormalizeData(data) #The user should not do this, we do it here to be sure!

    SizeBinned = floor.(Int64,size(data))
    BinnedData = zeros(T, SizeBinned...)

    Dimension = length(size(data)) #Get the dimension of the data, either 1D or 2D!
    for idx in CartesianIndices(BinnedData) 
        BoxIndex = floor.(Int64, (idx .- 1 ) / l ) .+ 1
        BinnedData[BoxIndex] += abs2.(data[idx])
    end
    return BinnedData
end


function ComputePartitionFunction(BinnedData::AbstractMatrix{T}, q::Number{V}) where {T<:Number,V<:Real}
    return sum(BinnedData .^ q)
end

function ComputePartitionFunction(data::AbstractMatrix{T}, qs::Vector{V}, l::Number) where {T<:Number,V<:Real}
    Zqs = zeros(T, length(qs))
    data = RenormalizeData(data) #The user should not do this, we do it here to be sure!
    BinnedData = BinData(data, l)
    for (i, q) in enumerate(qs)
        Zqs[i] = ComputePartitionFunction(BinnedData, q)
    end
    return Zqs
end

function ComputePartitionFunction(data::AbstractMatrix{T}, qs::Vector{T}, ls::Vector{V}) where {T<:Number,V<:Integer}
    data = RenormalizeData(data) #The user should not do this, we do it here to be sure!
    Zqs = zeros(T, length(qs), length(ls))
    for (j, l) in enumerate(ls)
        BinnedData = BinData(data, l)
        for (i, q) in enumerate(qs)
            Zqs[i, j] = ComputePartitionFunction(BinnedData, q)
        end
    end
    return Zqs
end

function ComputeMu(BinnedData::AbstractMatrix{T}, q::Number{V}) where {T<:Number,V<:Real} #Option if you don't have the partition function already
    return (BinnedData .^ q) ./ ComputePartitionFunction(BinnedData, q)
end

function ComputeMu(BinnedData::AbstractMatrix{T}, q::Number{V}, Zq::Number{T}) where {T<:Number,V<:Real} #option if you already have the partition function
    return (BinnedData .^ q) ./ Zq
end

function ObtainEntropy(μs::Vector{V}) where V<:Real
    return -sum(μs .* log.(μs) )
end


function ObtainZPrime(BinnedData::Vector{V}, q::Number{T}) where {T<:Real,V<:Real}
    return sum((BinnedData.^q) .* log.(BinnedData) )#.* (Ps .> 0))
end

function ObtainQs(qmin::Number{T}, qmax::Number{T}, num_q::Number{V}) where {T<:Real,V<:Integer}
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