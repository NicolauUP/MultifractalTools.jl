using MultifractalTools
using Test
using Aqua
using JET

@testset "MultifractalTools.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(MultifractalTools)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(MultifractalTools; target_defined_modules = true)
    end


    @testset "Core Logic & Monofractal Test" begin
    data = ones(Float64, 128, 128)
    qs = collect(LinRange(-2.0,2.0,11))
    scaling_data = compute_scaling_quantities(data, qs)
    
    n_scales = length(scaling_data.ls)
    @test n_scales > 0
    @test all(scaling_data.ls .<= scaling_data.ls[end])

    λ1 = 4
    λ2 = n_scales - 2
    spectrum_data = compute_spectrum(scaling_data, qs, λ1, λ2)

    @test length(spectrum_data.αs) == length(qs)
    f_max = maximum(spectrum_data.fs) #Must be ~ 2
    α_at_max = spectrum_data.αs[argmax(spectrum_data.fs)] #Must be ~ 2


    @test isapprox(f_max, 2.0; atol=0.01)
    @test isapprox(α_at_max, 2.0; atol=0.01)

    end
end