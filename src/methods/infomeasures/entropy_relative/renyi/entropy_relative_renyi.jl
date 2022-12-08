

# function entropy_relative(e::Renyi, est::Tsallis98, p::Probabilities, q::Probabilities)
#     q = e.q
#     return 1 / (q - 1) * (1 - sum(pᵢ^q / qᵢ^(q - 1) for (pᵢ, qᵢ) in zip(p, q)))
# end

# Analytical Renyi divergence for 1D normals:
# https://ieeexplore.ieee.org/abstract/document/6832827?casa_token=ZhfFH5_G6XgAAAAA:RzQMg0Zjn-CwtOWw4N-jeum3bWzP7tRioSSFAb76fZX58JmXDBW7mSqjbxvr73NDa9fplUSIGw
