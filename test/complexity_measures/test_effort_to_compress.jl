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
    end
end