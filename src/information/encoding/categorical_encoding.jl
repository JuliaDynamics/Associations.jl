import ComplexityMeasures: Encoding, encode, decode
export CategoricalEncoding, encode, decode

"""
    CategoricalEncoding <: Encoding
    CategoricalEncoding(x)

`CategoricalEncoding` is a generic encoding that encodes each `xᵢ ∈ unique(x)` to one of
the positive integers. Assumes that the elements of `x` are discrete.

The constructor requires the input data `x`, since the number of possible symbols
is `length(unique(x))`.
"""
struct CategoricalEncoding{T, I <: Integer} <: Encoding
    encode_dict::Dict{T, I}
    decode_dict::Dict{I, T}
    function CategoricalEncoding(x)

        # Ecode in order of first appearance, because `sort` doesn't work if we mix types,
        # e.g. `String` and `Int`.
        x_unique = unique(x)
        encode_dict = Dict{eltype(x), Int}()
        decode_dict = Dict{Int, eltype(x)}()
        for (i, xu) in enumerate(x_unique)
            encode_dict[xu] = i
            decode_dict[i] = xu
        end
        new{eltype(x), Int}(encode_dict, decode_dict)
    end
end
function CategoricalEncoding()
    throw(ArgumentError("`CategoricalEncoding` can't be initialized without input data."))
end

function encode(encoding::CategoricalEncoding, x)
    return encoding.encode_dict[x]
end

function decode(encoding::CategoricalEncoding, ω::Integer)
    return encoding.decode_dict[ω]
end
