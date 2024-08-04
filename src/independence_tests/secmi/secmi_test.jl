export SECMITest

"""
    SECMITest <: IndependenceTest
    SECMITest(; nshuffles = 19, surrogate = RandomShuffle(), rng = Random.default_rng())

A test for conditional independence based on the [`SECMI`](@ref) measure.
"""
Base.@kwdef struct SECMITest{E, S, I, RNG} <: IndependenceTest{E}
    est_or_measure::E = SECMI(; base = 2)
    surrogate::S = RandomShuffle()
    nshuffles::I = 19
    rng::RNG = Random.default_rng()
end

function independence(test::SECMITest, x, y, z)
    (; est_or_measure, surrogate, nshuffles, rng) = test
    sx = surrogenerator(x, surrogate, rng)
    secmi₀ = zeros(test.nshuffles)
    @show est_or_measure
    for k = 1:nshuffles
        @show typeof(sx()), typeof(y), typeof(z)
        secmi₀[k] = association(est_or_measure, sx(), y, z)
    end
    μ̂ = 1/nshuffles * sum(secmi₀)
    σ̂ = 1/(nshuffles - 1) * sum((sₖ - μ̂)^2 for sₖ in secmi₀)
    emp_cdf::Function = ecdf(secmi₀) 
    𝒩 = Normal(μ̂, σ̂)
    # degrees of freedom for Chi squared distribution estimated as the mean of the `secmiₖ`.
    # (page 18 in Kubkowski et al)
    𝒳 = Chisq(μ̂)
    D_N, D_chi2 = sup_values(emp_cdf, 𝒩, 𝒳, secmi₀)
    @show D_N, D_chi2
end

function sup_values(emp_cdf::Function, 𝒩, 𝒳, secmi₀)
    empirical_cdf_values = [emp_cdf(secmi₀, s) for s in secmi₀]
    normal_cdf_values = cdf.(𝒩, secmi₀)
    chi2_cdf_values = cdf.(𝒳, secmi₀)
    D_N = maximum(abs.(empirical_cdf_values .- normal_cdf_values))
    D_chi2 = maximum(abs.(empirical_cdf_values .- chi2_cdf_values))
    return D_N, D_chi2
end
