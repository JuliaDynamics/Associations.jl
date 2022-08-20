using Test

@testset "Effort to compress (ETC)" begin 
    @testset "ETC univariate" begin 
        alg = EffortToCompress(normalize = false)

        @testset "Paper examples" begin
            x1 = [1, 2, 1, 2, 1, 1, 1, 2]
            x2 = [0, 1, 0, 1, 0, 1, 0, 1]
            x3 = [0, 1, 0, 0, 1, 1, 1, 0]
            @test compression_complexity(x1, alg) == 5.0
            @test compression_complexity(x2, alg) == 1.0
            @test compression_complexity(x3, alg) == 6.0
        end

        @testset "MATLAB toolbox comparison" begin
            x4 = [0, 3, 2, 0, 1, 1, 1, 2, 2, 2, 1]
            x5 = [1, 0, 1, 1, 2, 1, 1, 1, 0, 2, 0]
            x6 = [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1]
            @test compression_complexity(x4, alg) == 10.0
            @test compression_complexity(x5, alg) == 8.0
            @test compression_complexity(x6, alg) == 6.0
        end
       
        @test compression_complexity([1], alg) == 0.0
        x3 = [0, 1, 0, 0, 1, 1, 1, 0]
        alg_norm = EffortToCompress(normalize = true)
        @test compression_complexity(x3, alg_norm) == 6.0 / (length(x3) - 1)

        # Sliding windows
        ts = rand(0:3, 50)
        alg = EffortToCompressSlidingWindow(window_size = 10, step = 2)
        windows = get_windows(ts, alg.window_size, alg.step)
        compression_complexity(ts, alg)
        @test length(compression_complexity(ts, alg)) == length(windows)
        @test all(compression_complexity(ts, alg) .>= 0)
    end

    @testset "ETC multivariate" begin
        x1 = [0, 0, 1, 0]
        x2 = [1, 0, 1, 0]
        x3 = [1, 1, 1, 1]
        x4 = [0, 1, 1, 1]

        @test symbol_sequence(Dataset(x1, x2, x3, x4), 2) == [6, 3, 15, 3]
        @test symbol_sequence(Dataset(x1, x3, x4), 2) == [2, 3, 7, 3]
        @test compression_complexity(Dataset(x1, x2, x3, x4), EffortToCompress(alphabet_size = 2)) == 3.0

        # Sliding window
        x1 = rand(0:1, 100)
        x2 = rand(0:1, 100)
        x3 = rand(0:1, 100)
        x4 = rand(0:1, 100)
        data = Dataset(x1, x2, x3, x4)
        window_size = 10
        step = 2
        alg = EffortToCompressSlidingWindow(
            alphabet_size = 2, 
            window_size = window_size, 
            step = step)
        
        windows = get_windows(data, window_size, step)
        res = compression_complexity(data, alg)
        @test length(res) == length(windows)
        @test all(res .>= 0)

        ###############################################
        # Multivariate with different dimensionalities
        ###############################################
        # Variables in the first dataset have an alphabet size of 2
        x1, x2 = rand(0:1, 1000), rand(0:1, 1000)

        # Variables in the second dataset have an alphabet size of 4
        y1, y2, y3 = rand(0:3, 1000), rand(0:3, 1000), rand(0:3, 1000)

        d1 = Dataset(x1, x2)
        d2 = Dataset(y1, y2, y3)
        alg = EffortToCompress(normalize = true)
        @test 0.0 .<= compression_complexity(d1, d2, alg, 2, 4) .<= 1.0

        windows = get_windows(data, 50, 10)
        alg = EffortToCompressSlidingWindow(normalize = true, window_size = 50, step = 10)
        sw = compression_complexity(data, alg)
        @test length(sw) == length(windows)
        @test all(0.0 .<= sw .<= 1.0)
    end

    @testset "ETC joint" begin
        # Test case in the "Joint ETC measure for a pair of time series: ETC(X,Y)" 
        # section in Kathpalia & Nagaraj (2013):
        alg = EffortToCompress(normalize = false)
        x = [1, 2, 1, 2, 1, 2]
        y = [1, 2, 1, 3, 1, 3]
        @test compression_complexity(x, y, alg) == 4.0

        # Constant sequences should not trigger compression
        x = [0, 0, 0, 0, 0]
        y = [1, 1, 1, 1, 1]
        compression_complexity(x, y, EffortToCompress()) == 0.0
       
        
        # Single-element sequences already have zero-entropy.
        x = [0]
        y = [1]
        compression_complexity(x, y, EffortToCompress()) == 0.0

        x = [1, 1, 0, 1, 1]
        y = [0, 0, 0, 0, 0]
        # Should reduce in three steps as
        # 11011 202 32 4
        # aaaaa bab cb d
        compression_complexity(x, y, EffortToCompress()) == 3.0

        # A good test of the correctness of the implementation is that the joint ETC(X, Y)
        # obeys the relationship stated in Kathpalia & Nagaraj (2013):
        # ETC(X, Y) <= ETC(X) + ETC(Y)
        # Here, we test if that is the case for our implementation by running 100000 test
        # cases on random pairs of length-50 sequences, whose number of symbols 
        # for each sequence is allowed to vary between 2 and 6 between iterations.
        inequality_tests = Vector{Bool}(undef, 0)
        alg = EffortToCompress(normalize = false)

        for i = 1:100000
            x = rand(0:rand(1:5), 50)
            y = rand(0:rand(1:5), 50)
            ETCxy = round(Int, compression_complexity(x, y, alg))
            ETCx = round(Int, compression_complexity(x, alg))
            ETCy = round(Int, compression_complexity(y, alg))
            res = ETCxy <= ETCx + ETCy
            push!(inequality_tests, res)
        end
        N_not_true = count(inequality_tests .!= true)
        N_true = count(inequality_tests .== true)
        @test all(inequality_tests .== true)

        @testset "ETC joint sliding window" begin
            x1 = rand(0:1, 500)
            x2 = rand(0:1, 500)

            window_size = 10
            step = 2
            alg = EffortToCompressSlidingWindow(
                window_size = window_size, 
                step = step)

            windows = get_windows(x1, window_size, step)
            res = compression_complexity(x1, x2, alg)
            @test length(res) == length(windows)
            @test all(res .>= 0)
        end

        @testset "ETC joint normalization" begin
            results = zeros(Float64, 10000)
            for i = 1:10000
                x, y = rand(0:3, 100), rand(0:5, 100)
                results[i] = compression_complexity(x, y, EffortToCompress(normalize = true))
            end
            @test all(0.0 .<= results .<= 1.0)
        end
    end
end