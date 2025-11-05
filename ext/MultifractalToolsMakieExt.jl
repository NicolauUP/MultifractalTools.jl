module MultifractalToolsMakieExt

# 1. Import your main package's functions
using MultifractalTools

# 2. Import the "heavy" dependencies this extension needs
using GLMakie

import LsqFit: curve_fit
import Statistics: mean

"""
    plot_to_fit(ScalingQuantities::NamedTuple, qs::AbstractVector{<:Real}) -> Figure

Launches an interactive Makie window to help visually determine the optimal
scaling region for the linear fits.

This plot shows \$\\log(Z(q, \\lambda))\$ vs \$\\log(\\lambda)\$. It provides sliders for:
- `q index`: To select which `q` value to inspect.
- `λ1 index (start)`: To set the starting point of the fit (green shaded region).
- `λ2 index (end)`: To set the ending point of the fit.

The plot automatically updates with a linear fit (dashed red line) and displays
the \$R^2\$ value for the selected region in the title, allowing for quick
assessment of the fit quality.

# Arguments
- `ScalingQuantities::NamedTuple`: The `NamedTuple` output from
  `compute_scaling_quantities`.
- `qs::AbstractVector`: The vector of `q` values used to generate `ScalingQuantities`.

# Returns
- `GLMakie.Figure`: The interactive Makie figure.

# Examples
```jldoctest
julia> # Create mock data for the test
julia> ls = [1, 2, 4, 8, 16];
julia> qs = [-5.0, 0.0, 5.0];
julia> Zqs = rand(5, 3) .+ 1.0;
julia> Sqs, ZPrimes = (rand(5, 3), rand(5, 3));
julia> scaling_data = (ls=ls, Zqs=Zqs, Sqs=Sqs, ZPrimes=ZPrimes);

julia> # The function is in the extension, so we must load GLMakie
julia> using GLMakie;

julia> fig = MultifractalTools.plot_to_fit(scaling_data, qs);
    
julia> isa(fig, Figure)
true
"""
function MultifractalTools.plot_to_fit(ScalingQuantities::NamedTuple, qs::AbstractVector{<:Real})

    # --- 1. Extract Data ---
    log_λs = log.(ScalingQuantities.ls ./ maximum(ScalingQuantities.ls)) 
    Zqs_matrix = ScalingQuantities.Zqs
    n_qs = length(qs)
    n_scales = length(ScalingQuantities.ls)

    # --- 2. Create the Figure and Slider ---
    fig = Figure(size = (700, 600))
    ax = Axis(fig[1, 1], 
        xlabel = L"\log(\lambda)", 
        ylabel = L"\log(Z(q,\lambda))",
        title = "Scaling Function for q = $(round(qs[1], digits=2))" # Initial title
    )
    
    # Right slider (q) - VERTICAL
    sg_right = SliderGrid(
        fig[2, 1], # Position: Middle-right
        (label = "q index", range = 1:n_qs, startvalue = n_qs,color_active = :indianred,color_active_dimmed=(:indianred, 0.5)),
        (label = "λ1 index (start)", range = 2:n_scales-1, startvalue = 2),
        (label = "λ2 index (end)", range = 2:n_scales-1, startvalue = n_scales-1)
    )
    # --- 3. Create the Observables ---
    q_idx_obs = sg_right.sliders[1].value 
    λ1_idx_obs = sg_right.sliders[2].value
    λ2_idx_obs = sg_right.sliders[3].value




    initial_log_Zq = log.(max.(Zqs_matrix[:, 1], 1e-100))
    log_Zq_obs = Observable([Point2f(l, z) for (l, z) in zip(log_λs[2:end-1], initial_log_Zq[2:end-1])])
    fit_line_obs = Observable(Point2f[])

# Observables for the vertical lines
    λ1_x_pos = @lift(log_λs[$λ1_idx_obs])
    λ2_x_pos = @lift(log_λs[$λ2_idx_obs])


    vspan!(ax, λ1_x_pos, λ2_x_pos, color = (:green, 0.15))



    # --- 4. Plot Everything ---
    # Scatter plot data is linked to the Observable
    scatter!(ax, log_Zq_obs, color=:blue, markersize=10, label="Raw Data") 
    lines!(ax, fit_line_obs, color=(:red, 0.5), linewidth=4, linestyle=:dash, label="Fit")



    # --- 5. THE LISTENER (Simplified) ---
    # This block runs ONLY when the q slider changes
    onany(q_idx_obs, λ1_idx_obs, λ2_idx_obs) do q_i, λ1, λ2
        
        # A. Calculate new scatter data
        new_Zq_data = log.(max.(Zqs_matrix[:, q_i], 1e-100))
        
        # B. PUSH the new data to the plot's observable
        log_Zq_obs[] = [Point2f(l, z) for (l, z) in zip(log_λs[2:end-1], new_Zq_data[2:end-1])]
        
        # C. Update the title
        ax.title = "Scaling Function for q = $(round(qs[q_i], digits=2))"

        ax.limits=(log_λs[2]-0.1, log_λs[end-1]+0.1, minimum(new_Zq_data[2:end-1])-0.1, maximum(new_Zq_data[2:end-1])+0.1)


        # C. Slice, Fit, and Calculate R²
        x_fit = log_λs[λ1:λ2]
        y_fit = new_Zq_data[λ1:λ2]

        fit_result = curve_fit(MultifractalTools.power_law_model, x_fit, y_fit, [1.0, 1.0])
        τ_q, intercept = fit_result.param

        r_squared = 1.0 - sum((y_fit .- MultifractalTools.power_law_model(x_fit, fit_result.param)).^2) / sum((y_fit .- mean(y_fit)).^2)


        fit_line_data = [Point2f(l, τ_q * l + intercept) for l in log_λs[λ1:λ2]]
        fit_line_obs[] = fit_line_data
        # D. Update the fit line in the plot
        ax.title = "R² = $(round(r_squared, digits=4))"


    end

    # --- 6. Trigger the listener once ---
    q_idx_obs[] = 1 
    
    return fig
end


"""
    plot_spectrum(SingularitySpectrumData::NamedTuple; which::Symbol = :Spectrum) -> Figure

Plots the calculated multifractal spectrum.

This function generates a plot of the singularity spectrum \$f(\\alpha)\$ vs \$\\alpha\$,
the mass scaling exponent \$\\tau(q)\$ vs \$q\$, or both, based on the `which` keyword.

# Arguments
- `SingularitySpectrumData::NamedTuple`: The `NamedTuple` output from `compute_spectrum`.

# Keyword Arguments
- `which::Symbol`: Controls which plot to generate.
    - `:Spectrum` (default): Plots the singularity spectrum, \$f(\\alpha)\$ vs \$\\alpha\$.
    - `:Tau`: Plots the multifractal exponent, \$\\tau(q)\$ vs \$q\$.
    - `:Both`: Generates a two-panel figure with both plots.

# Returns
- `GLMakie.Figure`: The Makie figure object containing the plot(s).

# Examples
```jldoctest
julia> # Create mock data for the test
julia> qs = -5:5.0;
julia> τqs = (qs .^ 2) ./ 10; # Mock tau(q)
julia> αs = qs ./ 5 .+ 0.5;   # Mock alpha(q)
julia> fs = 1.0 .- (αs .- 0.5).^2; # Mock f(alpha)
julia> spec_data = (qs=qs, τqs=τqs, αs=αs, fs=fs);

julia> # The function is in the extension, so we must load GLMakie
julia> using GLMakie;

julia> fig = MultifractalTools.plot_spectrum(spec_data, which=:Both);

julia> isa(fig, Figure)
true
"""
function MultifractalTools.plot_spectrum(SingularitySpectrumData::NamedTuple; which::Symbol = :Spectrum, SpectrumLimits = (0.0,3.0, 0.0, 3.0), TauLimits=(0.0,2.0, -2.0,2.0))


    plot_theme = theme_latexfonts()

    local fig

    if which == :Spectrum
        with_theme(plot_theme) do 
        fig = Figure(size = (800,800))
        ax = Axis(fig[1,1], xlabel = L"\alpha", ylabel = L"f(\alpha)", limits=SpectrumLimits)
        scatterlines!(ax, SingularitySpectrumData.αs, SingularitySpectrumData.fs, marker = :circle, markersize=12)

        end

    elseif which == :Tau
        with_theme(plot_theme) do 
        fig = Figure(size = (800,800))
        ax = Axis(fig[1,1], xlabel = L"q", ylabel = L"\tau(q)",limits=TauLimits)
        scatterlines!(ax, SingularitySpectrumData.qs, SingularitySpectrumData.τqs, marker = :circle, markersize=12)
        end

    elseif which == :Both
        with_theme(plot_theme) do 
        fig = Figure(size = (1600,800))
        ax1 = Axis(fig[1,1], xlabel = L"\alpha", ylabel = L"f(\alpha)", limits=SpectrumLimits)
        scatterlines!(ax1, SingularitySpectrumData.αs, SingularitySpectrumData.fs, marker = :circle, markersize=12)

        ax2 = Axis(fig[1,2], xlabel = L"q", ylabel = L"\tau(q)", limits=TauLimits)
        scatterlines!(ax2, SingularitySpectrumData.qs, SingularitySpectrumData.τqs, marker = :rect, markersize=12)


        end

    else
        error("which_type must be :Spectrum, :Tau or :Both")
    end
return fig
end


end