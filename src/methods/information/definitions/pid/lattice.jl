# ----------------------------------------------------------------------------------------
# The redundancy lattice of Williams & Beer (2010), for internal use.
#
# A node is an antichain α = {A₁, …, Aₖ} of nonempty subsets of the sources, none of which
# contains another. Nodes are ordered by
#     α ≼ β  ⟺  for every B ∈ β there is an A ∈ α with A ⊆ B.
# Given a redundancy function I∩(α; T) on the nodes, the partial information atoms Π are
# its Möbius inversion: I∩(α) = Σ_{β ≼ α} Π(β).
#
# Source subsets are stored as bitmasks: source i is bit i - 1.
# ----------------------------------------------------------------------------------------

const MAX_PID_SOURCES = 5

struct RedundancyLattice
    nsources::Int
    # Each node is a sorted vector of source-subset bitmasks.
    nodes::Vector{Vector{UInt}}
    index::Dict{Vector{UInt}, Int}
    # below[i, j] is true if node i is strictly below node j.
    below::BitMatrix
    # Node indices in an order where every node comes after all the nodes below it.
    order::Vector{Int}
end

is_submask(a::UInt, b::UInt) = (a & b) == a
mask_to_indices(mask::UInt) = [i for i in 1:(8 * sizeof(UInt)) if (mask >> (i - 1)) & 1 == 1]
indices_to_mask(idxs) = reduce(|, (UInt(1) << (i - 1) for i in idxs); init = UInt(0))

# All antichains of nonempty subsets of n sources. Subsets are added in increasing mask
# order, so each antichain is generated exactly once.
function lattice_antichains(n::Int)
    subsets = UInt.(1:(2^n - 1))
    result = Vector{Vector{UInt}}()
    current = UInt[]
    function extend!(start)
        for k in start:length(subsets)
            s = subsets[k]
            if all(c -> !is_submask(c, s) && !is_submask(s, c), current)
                push!(current, s)
                push!(result, sort(current))
                extend!(k + 1)
                pop!(current)
            end
        end
    end
    extend!(1)
    return result
end

# α ≼ β in the redundancy lattice.
lattice_precedes(α, β) = all(B -> any(A -> is_submask(A, B), α), β)

function RedundancyLattice(n::Int)
    1 ≤ n ≤ MAX_PID_SOURCES || throw(ArgumentError(
        "The redundancy lattice is implemented for 1 to $MAX_PID_SOURCES sources; got $n."))
    nodes = lattice_antichains(n)
    m = length(nodes)
    below = falses(m, m)
    for j in 1:m, i in 1:m
        if i != j && lattice_precedes(nodes[i], nodes[j])
            below[i, j] = true
        end
    end
    # If α ≺ β, everything below α is also below β, so α has strictly fewer nodes below
    # it. Sorting by that count therefore gives a valid order for the Möbius inversion.
    order = sortperm(vec(sum(below; dims = 1)))
    index = Dict(node => k for (k, node) in enumerate(nodes))
    return RedundancyLattice(n, nodes, index, below, order)
end

const REDUNDANCY_LATTICES = Dict{Int, RedundancyLattice}()
const REDUNDANCY_LATTICES_LOCK = ReentrantLock()

# The lattice for n sources, built once and cached (the 5-source lattice has 7579 nodes).
function redundancy_lattice(n::Int)
    lock(REDUNDANCY_LATTICES_LOCK) do
        get!(() -> RedundancyLattice(n), REDUNDANCY_LATTICES, n)
    end
end

# Partial information atoms from the redundancy values of every node.
function mobius_inversion(lattice::RedundancyLattice, redundancies::AbstractVector{T}) where T
    atoms = zeros(T, length(redundancies))
    m = length(redundancies)
    for j in lattice.order
        s = zero(T)
        for i in 1:m
            lattice.below[i, j] && (s += atoms[i])
        end
        atoms[j] = redundancies[j] - s
    end
    return atoms
end

# Convert a user-given antichain, e.g. `[[1], [2, 3]]`, to a node index.
function lattice_node_index(lattice::RedundancyLattice, α)
    n = lattice.nsources
    masks = UInt[]
    for A in α
        idxs = collect(Int, A)
        isempty(idxs) && throw(ArgumentError("Each element of an atom must be a nonempty set of sources."))
        all(i -> 1 ≤ i ≤ n, idxs) || throw(ArgumentError(
            "Source indices must be between 1 and $n; got $idxs."))
        push!(masks, indices_to_mask(idxs))
    end
    node = sort(unique(masks))
    haskey(lattice.index, node) || throw(ArgumentError(
        "$(collect(α)) is not an antichain: no element may contain another."))
    return lattice.index[node]
end

# Literature notation: {1}{23} for the antichain {{1}, {2, 3}}.
lattice_node_string(node) = join("{" * join(mask_to_indices(A)) * "}" for A in node)
lattice_node_vectors(node) = [mask_to_indices(A) for A in node]
