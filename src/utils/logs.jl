
# Just use ComplexityMeasures.convert_logunit when it is released.
"""
    _convert_logunit(h_a::Real, , to) → h_b

Convert a number `h_a` computed with logarithms to base `a` to an entropy `h_b` computed
with logarithms to base `b`. This can be used to convert the "unit" of e.g. an entropy.
"""
function _convert_logunit(h::Real, base_from, base_to)
    h / log(base_from, base_to)
end
