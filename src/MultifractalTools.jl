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
    return BinData
end


function ComputePartitionFunction(BinnedData::AbstractMatrix{T}, q::Number{V}) where {T<:Number,V<:Real}
    return sum(BinnedData .^ q)
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

end