
@testset "PredictiveAsymmetry" begin 
    x, y, z = rand(100), rand(100), rand(100)
    ηs = 1:3

    @testset "Binning based" begin

        @testset "VisitationFrequency" begin 
            # Define prediction lags and estimation method
            method = VisitationFrequency(RectangularBinning(3))

            # 𝔸(x → y) and  𝔸(x → y | z)
            𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
            𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
            𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
            𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

            @test 𝔸reg isa Vector{<:Real}
            @test 𝔸cond isa Vector{<:Real}
            @test 𝒜reg isa Vector{<:Real}
            @test 𝒜cond isa Vector{<:Real}
        end

        @testset "TransferOperator" begin 
            # Define prediction lags and estimation method
            method = TransferOperator(RectangularBinning(2))

            # 𝔸(x → y) and  𝔸(x → y | z)
            𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
            𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
            𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
            𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

            @test 𝔸reg isa Vector{<:Real}
            @test 𝔸cond isa Vector{<:Real}
            @test 𝒜reg isa Vector{<:Real}
            @test 𝒜cond isa Vector{<:Real}
        end
    end

    @testset "Nearest neighbor based" begin 
        @testset "Kraskov" begin 
            # Define prediction lags and estimation method
            method = Kraskov()

            # 𝔸(x → y) and  𝔸(x → y | z)
            𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
            𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
            𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
            𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

            @test 𝔸reg isa Vector{<:Real}
            @test 𝔸cond isa Vector{<:Real}
            @test 𝒜reg isa Vector{<:Real}
            @test 𝒜cond isa Vector{<:Real}
        end

        @testset "KozachenkoLeonenko" begin 
            # Define prediction lags and estimation method
            method = KozachenkoLeonenko()

            # 𝔸(x → y) and  𝔸(x → y | z)
            𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
            𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
            𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
            𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

            @test 𝔸reg isa Vector{<:Real}
            @test 𝔸cond isa Vector{<:Real}
            @test 𝒜reg isa Vector{<:Real}
            @test 𝒜cond isa Vector{<:Real}
        end
    end

    @testset "Kernel density based" begin 
        @testset "NaiveKernel" begin 
            # Define prediction lags and estimation method
            method = NaiveKernel(0.2)

            # 𝔸(x → y) and  𝔸(x → y | z)
            𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
            𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
            𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
            𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

            @test 𝔸reg isa Vector{<:Real}
            @test 𝔸cond isa Vector{<:Real}
            @test 𝒜reg isa Vector{<:Real}
            @test 𝒜cond isa Vector{<:Real}
        end
    end


    @testset "Hilbert" begin
        # Define prediction lags and estimation method
        method = Hilbert(VisitationFrequency(RectangularBinning(3)))

        # 𝔸(x → y) and  𝔸(x → y | z)
        𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
        𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
        𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
        𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

        @test 𝔸reg isa Vector{<:Real}
        @test 𝔸cond isa Vector{<:Real}
        @test 𝒜reg isa Vector{<:Real}
        @test 𝒜cond isa Vector{<:Real}
    end

    @testset "Symbolic" begin 
        @testset "SymbolicPermutation" begin 
            # Define prediction lags and estimation method
            method = SymbolicPermutation(m = 3)

            # 𝔸(x → y) and  𝔸(x → y | z)
            𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
            𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
            𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
            𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

            @test 𝔸reg isa Vector{<:Real}
            @test 𝔸cond isa Vector{<:Real}
            @test 𝒜reg isa Vector{<:Real}
            @test 𝒜cond isa Vector{<:Real}
        end
    end

end