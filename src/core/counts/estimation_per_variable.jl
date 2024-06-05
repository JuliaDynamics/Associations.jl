using ComplexityMeasures: encode

# If only one encoding is given, apply same encoding to all points
function counts(encoding::VariableEncoding{1}, x::Vararg{ArrayOrStateSpaceSet, N}) where N
    e = first(encoding.encodings)
    x̂ = ([encode(e, pt) for pt in xₖ] for xₖ in x)
    return counts(x̂...)
end
