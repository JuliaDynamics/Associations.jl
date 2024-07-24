# [Encoding](@id encoding_api_and_tutorial)

To estimate [`Probabilities`](@ref), which are the input to [`MultivariateInformationMeasure`](@ref)s, 
we encode "encode" input data into an intermediate representation indexed by the positive integers. 
This intermediate representation is called an "encoding".

We here re-export relevant types and functions [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl) that perform this type of coarse-graining.

These encoding schemes are used as input to [`CodifyPoints`](@ref).

```@docs
Encoding
GaussianCDFEncoding
OrdinalPatternEncoding
RelativeMeanEncoding
RelativeFirstDifferenceEncoding
UniqueElementsEncoding
RectangularBinEncoding
CombinationEncoding
encode
decode
```