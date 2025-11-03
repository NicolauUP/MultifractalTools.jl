module MultifractalTools

using LinearAlgebra
using LsqFit
using Statistics
using GLMakie



export obtain_qs, compute_scaling_quantities, compute_spectrum, plot_spectrum, get_divisors, find_best_scaling_size


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

function plot_to_fit(ScalingQuantities::NamedTuple, qs::AbstractVector{<:Real})

    λs = log.(ScalingQuantities.ls ./ maximum(ScalingQuantities.ls))

    Zqs = ScalingQuantities.Zqs

    n_qs = length(qs)
    n_λs = length(λs)

    # Create the Figure and Sliders

    fig = Figure(size = (900, 700))
    ax = Axis(fig[2,1], xlabel = L"\log(\lambda)", ylabel = L"\log(Z(q,\lambda))")

    sg_top = SliderGrid(
        fig[1,1],
        (label = "Select starting index",range = 1:n_scales, startvalue = 1),
        (label = "Select ending index", range = 1:n_scales, startvalue = n_scales)
    )

    sg_right = SliderGrid(
        fig[2,2],
        (label = "Select q index", range = 1:n_qs, startvalue = 1)
    )


    #Create the Observables!

    q_idx = sg_right.sliders[1].values
    λ_start_idx = sg_top.sliders[1].values
    λ_end_idx = sg_top.sliders[2].values

    
    log_Zq_obs = @lift(log.(ScalingQuantities.Zqs[1:end, $q_idx]))

    fit_results = @lift begin
        q_i = $q_idx
        λ1 = $λ_start_idx
        λ2 = $λ_end_idx
         
        if λ1 >= λ2
            (title = "Invalid λ range", line = Point2f[])
        else

            #. 1. Get data
            log_Zq_data = $log_Zq_obs


            # 2. Slice for fitting
            x_fit = log_λs[λ1:λ2]
            y_fit = log_Zq_data[λ1:λ2]

            # 3. Fit
            fit = curve_fit(power_law_model, x_fit, y_fit, [1.0, 1.0])
            τ_q = fit.param[1]
            intercept = fit.param[2]

            # 4. Get R^2
            residulas = y_fit  .- power_law_model(x_fit, fit.param)
            R2 = 1.0 - sum(residulas .^ 2) / sum((y_fit .- mean(y_fit)) .^ 2)

            # 5. Create line data for the plot!
            fit_line_data = [Point2f(l, τ_q*l + intercept) for l in log_λs[λ1:λ2]]


            # 6. Create title
            title_str = "R^2 = $(round(R2, digits=4)), q = $(qs[q_i])"

            (title = title_str, line = fit_line_data)
        end
    end

    # ---- 5. Plot -----

    scatter!(ax, log_ls, log_Zq_obs, marker = :circle, markersize=10)

    lines!(ax, @lift(fit_results.line), color = :black, linewidth = 2,linestyle=:dash)

    ax.title = @lift(fit_results.title)
    
    set_close_to!(sg_right.slides[1],1)

    return fig
end

function plot_spectrum(SingularitySpectrumData::NamedTuple; which_type::Symbol = :Spectrum)


    if which_type == :Spectrum
        with_theme(theme_latexfonts()) do 
        fig = Figure(size = (800,800))
        ax = Axis(fig[1,1], xlabel = L"\alpha", ylabel = L"f(\alpha)")
        scatterlines!(ax, SingularitySpectrumData.αs, SingularitySpectrumData.fs, marker = :circle, markersize=12)
        display(fig)
        end

    elseif which_type == :Tau
        with_theme(theme_latexfonts()) do 
        fig = Figure(size = (800,800))
        ax = Axis(fig[1,1], xlabel = L"q", ylabel = L"\tau(q)")
        scatterlines!(ax, SingularitySpectrumData.qs, SingularitySpectrumData.τqs, marker = :circle, markersize=12)
        display(fig)
        end

    elseif which_type == :Both
        with_theme(theme_latexfonts()) do 
        fig = Figure(size = (1600,800))
        ax1 = Axis(fig[1,1], xlabel = L"\alpha", ylabel = L"f(\alpha)")
        scatterlines!(ax1, SingularitySpectrumData.αs, SingularitySpectrumData.fs, marker = :circle, markersize=12)

        ax2 = Axis(fig[1,2], xlabel = L"q", ylabel = L"\tau(q)")
        scatterlines!(ax2, SingularitySpectrumData.qs, SingularitySpectrumData.τqs, marker = :rect, markersize=12)

        display(fig)
        end

    else
        error("which_type must be :Spectrum, :Tau or :Both")
    end

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