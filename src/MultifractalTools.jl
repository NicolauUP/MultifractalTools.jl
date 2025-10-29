module MultifractalTools

using LinearAlgebra
using LsqFit
using Statistics


function RenormalizeData(data::AbstractMatrix{T}) where T<:Number
    value = sum(abs2.(data))
    return data / sqrt(value)
end
