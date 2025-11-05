module MultifractalTools

using LinearAlgebra
using LsqFit
using Statistics



export obtain_qs, compute_scaling_quantities, compute_spectrum, plot_spectrum, plot_to_fit



"""
    renormalize_data(data::AbstractMatrix)

Normalizes the input `data` matrix according to its \$L^2\$-norm.

This ensures that the sum of the absolute squares of the elements is equal to 1.
Mathematically, this computes:
\$D' = D / \\sqrt{\\sum |D_{ij}|^2}\$

This is typically the first step in a multifractal analysis to treat the
data as a probability distribution (\$\\|\\psi\\|^2 = 1\$).

# Arguments
- `data::AbstractMatrix`: The input 2D data.

# Returns
- A new `AbstractMatrix` of the same size, with its elements scaled.

# Examples
```jldoctest
julia> A = [1.0 1.0; 1.0 1.0]
2×2 Matrix{Float64}:
 1.0  1.0
 1.0  1.0

julia> B = MultifractalTools.renormalize_data(A);

julia> sum(abs2.(B))
1.0
"""
function renormalize_data(data::AbstractMatrix{T}) where T<:Number
    value = sum(abs2.(data))
    return data / sqrt(value)
end

"""
    get_divisors(n::Integer) -> Vector{Int}

Finds all positive divisors of an integer `n` and returns them in a sorted vector.

The function computes divisors by iterating up to `sqrt(n)` for efficiency.
It correctly handles `n=0` (returning an empty list) and negative numbers
(by taking the absolute value).

# Arguments
- `n::Integer`: The number to find the divisors of.

# Returns
- A `Vector{Int}` containing all divisors of `n` in ascending order.

# Examples
```jldoctest
julia> MultifractalTools.get_divisors(100)
9-element Vector{Int64}:
   1
   2
   4
   5
  10
  20
  25
  50
 100

julia> MultifractalTools.get_divisors(17)
2-element Vector{Int64}:
  1
 17
"""
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


"""
    find_best_scaling_size(max_size::Integer, crop_ratio::Float64) -> NamedTuple

Finds an optimal square size for box-counting analysis near `max_size`.

This function searches for an integer `n` in the range 
`[floor(max_size * (1 - crop_ratio)), max_size]` that has the maximum number 
of divisors. This is useful for maximizing the number of available scaling 
sizes (`l` values) for the analysis, which improves the quality of the 
linear fits.

# Arguments
- `max_size::Integer`: The largest possible size (e.g., `minimum(size(data))`).
- `crop_ratio::Float64`: The fraction of the `max_size` to search within. 
  For example, `0.1` searches in the top 10% of the range.

# Returns
- A `NamedTuple` with two fields:
    - `size`: The integer `n` found to have the most divisors.
    - `divisors`: A `Vector{Int}` of the divisors of `n`.

# Examples
```jldoctest
# Search in the range [90, 100]
# 96 has 12 divisors: [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 96]
# 100 has 9 divisors: [1, 2, 4, 5, 10, 20, 25, 50, 100]
julia> result = MultifractalTools.find_best_scaling_size(100, 0.1);

julia> result.size
96

julia> result.divisors
12-element Vector{Int64}:
  1
  2
  3
  4
  6
  8
 12
 16
 24
 32
 48
 96
"""
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

"""
    bin_data(data::AbstractMatrix, l::Integer) -> Matrix

Performs spatial coarse-graining (binning) of the data.

This function divides the input `data` matrix into non-overlapping square 
boxes of size `l x l`. It calculates the "measure" (probability) in each 
box by summing the `abs2` of all elements within that box.

The `abs2` (e.g., `|ψ|^2`) is used, as the multifractal measure `μ` is 
typically a probability distribution.

# Arguments
- `data::AbstractMatrix`: The input 2D data, assumed to be renormalized 
  (though not strictly required).
- `l::Integer`: The side length of the square bins.

# Returns
- A new `Matrix` (named `BinnedData`) of size `floor.(size(data) ./ l)`, 
  where each element is the sum of `abs2(data[...])` within that bin.

# Examples
```jldoctest
julia> A = [1.0 1.0 2.0 2.0; 
            1.0 1.0 2.0 2.0; 
            3.0 3.0 4.0 4.0; 
            3.0 3.0 4.0 4.0]
4×4 Matrix{Float64}:
 1.0  1.0  2.0  2.0
 1.0  1.0  2.0  2.0
 3.0  3.0  4.0  4.0
 3.0  3.0  4.0  4.0

# Bin with l=2
# Box 1 (top-left): 1²+1²+1²+1² = 4
# Box 2 (top-right): 2²+2²+2²+2² = 16
# Box 3 (bottom-left): 3²+3²+3²+3² = 36
# Box 4 (bottom-right): 4²+4²+4²+4² = 64
julia> MultifractalTools.bin_data(A, 2)
2×2 Matrix{Float64}:
  4.0  16.0
 36.0  64.0
"""
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


"""
    compute_scaling_quantities(data, qs; ls=Integer[], crop_to_best_fit=true, crop_ratio=0.1)

Computes the core scaling quantities required for multifractal analysis.

This function calculates the partition function \$Z(q, l)\$, the entropy \$S(q, l)\$, 
and the average \$\\log(\\mu)\$ weighted by \$\\mu^q\$ (which we call \$Z'(q, l)\$) for a given 
`data` matrix, across a range of `q` values and box sizes `l`.

# Arguments
- `data::AbstractMatrix`: The input 2D data (e.g., a wavefunction or probability distribution).
- `qs::AbstractVector`: A vector of `q` values (the "moments") to analyze.

# Keyword Arguments
- `ls::AbstractVector{<:Integer}`: An optional vector of box sizes `l` to use. If left empty 
  (default), the function will automatically determine box sizes based on the `data` size 
  and other keyword arguments.
- `crop_to_best_fit::Bool`: If `true` (default) and `ls` is empty, the function will find 
  the largest "best" square size (with the most divisors) within `crop_ratio` of the 
  original data size. If `false`, it uses the full data size.
- `crop_ratio::Float64`: The fraction of the data size to search for a "best" size 
  (default: `0.1`). For a 1000x1000 matrix, it will search from 900-1000.

# Returns
- A `NamedTuple` containing:
    - `ls`: The vector of box sizes `l` that were used.
    - `Zqs`: An `(l, q)` matrix of the partition functions, \$Z(q, l) = \\sum_i \\mu_i^q\$.
    - `Sqs`: An `(l, q)` matrix of the information entropies, \$S(q, l) = \\sum_i p_i \\log(p_i)\$.
    - `ZPrimes`: An `(l, q)` matrix of the weighted log averages, \$Z'(q, l) = \\sum_i p_i \\log(\\mu_i)\$.

# Examples
```jldoctest
julia> data = rand(128, 128);
julia> qs = obtain_qs(-5, 5, 11);
julia> scaling_data = compute_scaling_quantities(data, qs; crop_to_best_fit=false);
julia> scaling_data.ls
8-element Vector{Int64}:
   1
   2
   4
   8
  16
  32
  64
 128

"""
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

"""
    power_law_model(x, p)

A simple linear model \$f(x) = p_1 x + p_2\$ used for fitting power-law data.

This function is intended to be passed to `LsqFit.curve_fit` to find scaling
exponents from log-transformed data. 

It represents a power law \$y = C x^\\alpha\$ that has been transformed into
log-log space: \$\\log(y) = \\alpha \\log(x) + \\log(C)\$.

In this model:
- `x` corresponds to the log-transformed independent variable (e.g., \$\\log(l)\$).
- `p[1]` is the parameter for the slope (the scaling exponent \$\\alpha\$).
- `p[2]` is the parameter for the y-intercept (the prefactor \$\\log(C)\$).

# Arguments
- `x`: The independent variable (e.g., `log.(ls)`).
- `p`: A 2-element vector of parameters `[slope, intercept]`.

# Returns
- The fitted value `p[1] .* x .+ p[2]`.

# Examples
```jldoctest
julia> using LsqFit

# Create data for a line y = 2x + 1
julia> x_data = [1.0, 2.0, 3.0, 4.0];
julia> y_data = [3.0, 5.0, 7.0, 9.0];
julia> p0 = [1.0, 0.0]; # Initial guess for [slope, intercept]

julia> fit = curve_fit(MultifractalTools.power_law_model, x_data, y_data, p0);
    
julia> fit.param
2-element Vector{Float64}:
 2.0
 1.0
""" 
function power_law_model(x, p) #Log Scale
    return p[1] .* x .+ p[2]
end

"""
    compute_spectrum(ScalingQuantities, qs, λ1, λ2) -> NamedTuple

Calculates the multifractal spectrum exponents (\$\tau(q)\$, \$\alpha(q)\$, \$f(\alpha)\$)
by fitting the results from `compute_scaling_quantities`.

This function performs a linear regression in log-log space to find the
scaling exponents. For example, the multifractal exponent \$\tau(q)\$ is found by
fitting the power law \$Z(q, \epsilon) \\sim \epsilon^{\tau(q)}\$, which becomes
linear in log space: \$\\log(Z) = \tau(q) \\log(\epsilon) + C\$.

The function fits this linear relationship using the `power_law_model`
for \$\tau(q)\$, \$\alpha(q)\$, and \$f(\alpha)\$ simultaneously.

# Arguments
- `ScalingQuantities::NamedTuple`: The `NamedTuple` output from
  `compute_scaling_quantities`.
- `qs::AbstractVector`: The vector of `q` values. This must be the same
  `qs` vector used to generate the `ScalingQuantities`.
- `λ1::Integer`: The **starting index** (not value) from `ScalingQuantities.ls`
  to use for the linear fit.
- `λ2::Integer`: The **ending index** from `ScalingQuantities.ls` to use
  for the linear fit.

# Returns
- A `NamedTuple` containing the calculated spectra:
    - `qs`: The `q` values.
    - `τqs`: The vector of multifractal exponents, \$\tau(q)\$.
    - `αs`: The vector of singularity spectrum \$\alpha(q)\$.
    - `fs`: The vector of singularity spectrum \$f(\alpha)\$.

# Examples
```jldoctest
julia> data = rand(128, 128);
julia> qs = MultifractalTools.obtain_qs(-5, 5, 3); # -5.0, 0.0, 5.0
julia> scaling_data = MultifractalTools.compute_scaling_quantities(data, qs);
julia> scaling_data.ls
8-element Vector{Int64}:
   1
   2
   4
   8
  16
  32
  64
 128

# We fit using the scaling range from index 2 (l=2) to index 6 (l=32)
julia> spectrum = compute_spectrum(scaling_data, qs, 2, 6);
julia> spectrum.τqs
3-element Vector{Float64}:
  1.9999999999999996
 -0.0
 -2.0
"""
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



"""
    obtain_qs(qmin::Number, qmax::Number, num_q::Integer) -> Vector{Float64}

Generates a vector of `num_q` `q`-values, linearly spaced between `qmin` and `qmax`.

# Arguments
- `qmin`: The minimum value for the `q` range.
- `qmax`: The maximum value for the `q` range.
- `num_q`: The number of `q` values to generate.

# Returns
- A `Vector` of `Float64` values.

# Examples
```jldoctest
julia> obtain_qs(-5, 5, 3)
3-element Vector{Float64}:
 -5.0
  0.0
  5.0
"""
function obtain_qs(qmin::Number, qmax::Number, num_q::Integer) 
    return collect(LinRange(qmin, qmax, num_q)) #This could be memory-problematic but for now it's ok
end


end #module