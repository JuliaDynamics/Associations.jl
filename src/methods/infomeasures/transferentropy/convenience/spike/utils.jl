using StaticArrays, DelayEmbeddings

"""
    EventIdentifier

An abstract type indicating some sort of event.
"""
abstract type EventIdentifier end

"""
    OneSpike <: EventIdentifier
    OneSpike()

Events are indicated by the presence of `1`s or `true`s.

## Example

```julia
# Events occur at indices 2, 4, and 6.
x = [0, 1, 0, 1, 0, 1]

# Events again occur at indices 2, 4, and 6.
x = [-2, 1, 5, 1, 2, 1]
```
"""
struct OneSpike <: EventIdentifier end

"""
    first_index(x::AbstractVector{T}, m::Int) where T

Find the first time index for which a continuous-time embedding with history
length `m` can be constructed.
"""
function first_index(x::AbstractVector{T}, event::EventIdentifier = OneSpike();
        m::Int = 2) where T

    spike = one(T)

    i::Int = 1
    ct::Int = 0
    while ct < m && i < length(y)
        if x[i] == spike
            ct += 1
        end
        i += 1
    end

    ct >= m || throw(ErrorException("Could not find m=$m spikes in `x`"))

    return i
end

function first_index(x::AbstractDataset{D, T}, event::EventIdentifier = OneSpike();
        m::Int = 2) where {D, T}

    spike = one(T)
    L = length(x)
    cts::MVector{D, Int} = zeros(MVector{D, Int})

    i::Int = 1
    while any(cts .< m) && i < L
        for d = 1:D
            if x[i][d] == spike
                cts[d] += 1
            end
        end
        i += 1
    end

    all(cts .>= m) ||
        throw(ErrorException("Couldn't find m=$m spikes in one or more columns of `x`"))
    return i
end

first_index(xs; kwargs...) = map(v -> first_index.(v; kwargs...), xs)

# using Test
# using DelayEmbeddings
# z = [1, 0, 1, 0, 0, 1, 1, 0, 1, 0]
# x = [0, 1, 0, 1, 0, 1, 0, 0, 0, 1]
# y = [1, 1, 0, 0, 1, 0, 1, 0, 0, 1]
# D = Dataset(x, y, z)
# @test first_index(z, m = 2) == 4
# @test first_index(z, m = 3) == 7
# @test first_index(z, m = 4) == 8
# @test first_index(D, m = 2) == 5
# @test first_index(D, m = 3) == 7
# @test_throws ErrorException first_index(D, m = 4)
