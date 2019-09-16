
# [TEVars](@id TEVars)

A `TEVars` instance maps variables of a [custom delay reconstruction](@ref custom_delay_reconstruction) to the correct marginals during [TE computation](@ref general_workflow_te). In practice, this translates to using the correct columns of the `Dataset` furnishing the delay reconstruction when computing the different marginals 
during transfer entropy computation.

## Map reconstruction variables to marginals using keywords

```@docs 
TEVars()
```

## Infer the mapping between reconstruction variables and marginals

Alternatively, the mapping is inferred from the order of the inputs. Here, `target_future` is the ``T_{f}`` component, `target_presentpast` is the ``T_{pp}`` component, and `source_presentpast` is the 
``S_{pp}`` component of the delay reconstruction. For conditional analyses, the additional 
component `conditioned_presentpast` (the ``C_{pp}`` component) is needed.

```@docs 
TEVars(::Vector{Int}, ::Vector{Int}, ::Vector{Int})
```

```@docs 
TEVars
```