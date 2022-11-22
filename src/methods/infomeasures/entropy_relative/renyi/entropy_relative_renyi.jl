

# function entropy_relative(e::Renyi, est::Tsallis98, p::Probabilities, q::Probabilities)
#     q = e.q
#     return 1 / (q - 1) * (1 - sum(pᵢ^q / qᵢ^(q - 1) for (pᵢ, qᵢ) in zip(p, q)))
# end
