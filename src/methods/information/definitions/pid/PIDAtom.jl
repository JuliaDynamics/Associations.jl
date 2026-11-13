export PIDAtom

"""
    PIDAtom <: MultivariateInformationMeasure
    PIDAtom(pid::PartialInformationDecomposition, atom)

A single atom of a [`PartialInformationDecomposition`](@ref), as a number. `atom` is one of

- `:synergy`: the atom ``\\{12 \\ldots n\\}`` (see [`synergy`](@ref));
- `:redundancy`: the atom ``\\{1\\}\\{2\\} \\ldots \\{n\\}`` (see [`redundancy`](@ref));
- `:unique => i`: the unique information of source `i` (see [`unique_information`](@ref));
- a lattice node given as a collection of source sets, e.g. `[[1], [2, 3]]`.

## Usage

- Use with [`association`](@ref) to compute the atom from input data:
    `association(est, s₁, s₂, …, target)` returns a number. The last input is the target.
- Use with [`SurrogateAssociationTest`](@ref) to test whether the atom is significantly
    larger than expected by chance. The test shuffles the **first** input, so put the
    source of interest first; any atom that involves source 1 (including synergy and
    redundancy) is zero when source 1 is independent of the other inputs.

Compatible estimators are those of the wrapped [`PartialInformationDecomposition`](@ref).

## Examples

```julia
using Associations
x, y = rand(0:1, 1000), rand(0:1, 1000)
t = x .⊻ y
pid = PartialInformationDecomposition(IMMI())
est = JointProbabilities(PIDAtom(pid, :synergy), CodifyVariables(UniqueElements()))
association(est, x, y, t) # ≈ 1 bit
independence(SurrogateAssociationTest(est), x, y, t)
```
"""
struct PIDAtom{P<:PartialInformationDecomposition, A} <: MultivariateInformationMeasure
    pid::P
    atom::A
    function PIDAtom(pid::P, atom::A) where {P<:PartialInformationDecomposition, A}
        if atom isa Symbol
            atom in (:synergy, :redundancy) || throw(ArgumentError(
                "Atom symbol must be `:synergy` or `:redundancy`; got :$atom."))
        elseif atom isa Pair
            (first(atom) == :unique && last(atom) isa Integer) || throw(ArgumentError(
                "A pair atom must be `:unique => i`; got $atom."))
        end
        return new{P, A}(pid, atom)
    end
end

min_inputs_vars(a::PIDAtom) = min_inputs_vars(a.pid)
max_inputs_vars(a::PIDAtom) = max_inputs_vars(a.pid)

extract_pid_atom(res::PIDResult, a::Symbol) = a == :synergy ? synergy(res) : redundancy(res)
extract_pid_atom(res::PIDResult, a::Pair) = unique_information(res, last(a))
extract_pid_atom(res::PIDResult, a) = atom(res, a)

# The same estimator with another definition.
with_definition(est::JointProbabilities, def) = JointProbabilities(def, est.discretization, est.pest)
with_definition(est::EntropyDecomposition{<:Any, <:DifferentialInfoEstimator}, def) =
    EntropyDecomposition(def, est.est)
with_definition(est::EntropyDecomposition{<:Any, <:DiscreteInfoEstimator}, def) =
    EntropyDecomposition(def, est.est, est.discretization, est.pest)
with_definition(est::MIDecomposition, def) = MIDecomposition(def, est.est)

function pid_atom_association(est, x...)
    a = est.definition
    res = association(with_definition(est, a.pid), x...)
    return extract_pid_atom(res, a.atom)
end

const PID_ATOM_ESTIMATORS = Union{
    JointProbabilities{<:PIDAtom},
    EntropyDecomposition{<:PIDAtom},
    MIDecomposition{<:PIDAtom},
}
association(est::PID_ATOM_ESTIMATORS, x...) = pid_atom_association(est, x...)
