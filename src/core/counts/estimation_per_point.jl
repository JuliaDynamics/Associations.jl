
# If multiple encodings are given, the number of encodings must match the number of
# input variables.
function counts(encoding::PointEncoding{N}, x::Vararg{ArrayOrStateSpaceSet, N}) where N
    x̂ = encode(encoding, x...)
    return counts(x̂...)
end


# If only one encoding is given, apply same encoding to all points
function counts(encoding::PointEncoding{1}, x::Vararg{ArrayOrStateSpaceSet, N}) where N
    e = first(encoding.encodings)
    x̂ = ([encode(e, pt) for pt in xₖ] for xₖ in x)
    return counts(x̂...)
end
