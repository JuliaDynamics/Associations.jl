
x = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0]
y = [0, 0, 0, 1, 0, 0, 1, 0, 0, 1]
export lean

"""
    lean_ex(C, E, l)

Computes the leaning of the `l`-assignment ``\\{C, E\\} = \\{X, Y\\}, \\{X_{t-l}, Y_t\\}``.
"""
function lean(X, Y, l = 1)
    @assert length(X) == length(Y)
    L = length(X)
    
    # Leave reasonably many for probability computations. `l` can't be too large.
    @assert l < L ÷ 2

    #P(y_t = 1 | x_{t-1} = 1) = n_ec/n_c

    # Number of time the state y_t = 1 AND x_{t-1} = 1 appears in {X, Y}
    n_ecs = Vector{Int}(undef, 0)
    n_cs = Vector{Int}(undef, 0)
    n_es = Vector{Int}(undef, 0)

    states = unique([X; Y])
    
    # Define an iterator over all penchants
    penchants = Iterators.product(states, states)
    @show states
    @show penchants |> collect
    

      

    for penchant in penchants
        @show penchant
        sx, sy = penchant[1], penchant[2]
            
        # Number of times the assumed cause has appeared
        n_c = 0

        # Number of times the assumed effect has appeared
        n_e = 0

        # Number of time the state y_t = a AND x_{t-1} = b appears in {X, Y},
        # where (a, b) is one of the length(states)^2 possible states.
        n_ec = 0
        for t = l+1:L
            effect = Y[t]
            cause = X[t-l]

            if effect == sx && cause == sy
                n_ec += 1
            end
            

            if effect == sx
                n_e += 1
            end

            if cause == sy
                n_c += 1
            end
        end
        push!(n_ecs, n_ec)
        push!(n_cs, n_c)
        push!(n_es, n_e)
    end

    @show Ps_cond = n_ecs ./ n_cs
    @show Ps_e = n_es ./ (L - l)
    @show Ps_c = n_cs ./ (L - l)
    
    # for t = l+1:L
    #     yt = Y[t]
    #     xtl = X[t-l]
    #     @show Y[t], X[t-l]
    #     if yt == 1 && xtl == 1
    #         n_ec += 1
    #     end

    #     if X[t-l] == 1
    #         n_c += 1
    #     end
    # end

    # @show n_ec
    # @show n_c

    # p = n_ec/n_c
end

