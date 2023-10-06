
# If multiple encodings are given, the number of encodings must match the number of
# input variables.
function counts(encoding::CodifyPoints{N}, x::Vararg{T, N}) where {T, N}
    x̂ = codify(encoding, x...)
    return counts(UniqueElements(), x̂...)
end


# If only one encoding is given, apply same encoding to all points
function counts(encoding::CodifyPoints{1}, x::Vararg{T, N}) where {T, N}
    e = first(encoding.encodings)
    x̂ = ([encode(e, pt) for pt in xₖ] for xₖ in x)
    return counts(UniqueElements(), x̂...)
end
