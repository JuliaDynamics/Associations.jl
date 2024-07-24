using DSP
using ComplexityMeasures: EntropyDefinition

export Amplitude, Phase, Hilbert

abstract type InstantaneousSignalProperty end

"""
    Amplitude <: InstantaneousSignalProperty

Indicates that the instantaneous amplitudes of a signal should be used. """
struct Amplitude <: InstantaneousSignalProperty end

"""
    Phase <: InstantaneousSignalProperty

Indicates that the instantaneous phases of a signal should be used. """
struct Phase <: InstantaneousSignalProperty end

# TODO: update to new syntax
"""
    Hilbert(est;
        source::InstantaneousSignalProperty = Phase(),
        target::InstantaneousSignalProperty = Phase(),
        cond::InstantaneousSignalProperty = Phase())
    ) <: TransferDifferentialEntropyEstimator

Compute transfer entropy on instantaneous phases/amplitudes of relevant signals, which are
obtained by first applying the Hilbert transform to each signal, then extracting the
phases/amplitudes of the resulting complex numbers [Palus2014](@cite). Original time series are
thus transformed to instantaneous phase/amplitude time series. Transfer
entropy is then estimated using the provided `est` on those phases/amplitudes (use e.g.
[`VisitationFrequency`](@ref), or [`OrdinalPatterns`](@ref)).

!!! info
    Details on estimation of the transfer entropy (conditional mutual information)
    following the phase/amplitude extraction step is not given in Palus (2014). Here,
    after instantaneous phases/amplitudes have been obtained, these are treated as regular
    time series, from which transfer entropy is then computed as usual.

See also: [`Phase`](@ref), [`Amplitude`](@ref).
"""
struct Hilbert{M} <: TransferEntropyEstimator{M}
    est::M # the estimator of the transfer entropy to be applied *after* transformation of the time series.
    source::InstantaneousSignalProperty
    target::InstantaneousSignalProperty
    cond::InstantaneousSignalProperty

    function Hilbert(est::M;
            source::InstantaneousSignalProperty = Phase(),
            target::InstantaneousSignalProperty = Phase(),
            cond::InstantaneousSignalProperty = Phase()) where M
        new{M}(est, source, target, cond)
    end
end

function association(est::Hilbert, source, target)
    hil_s = DSP.hilbert(source)
    hil_t = DSP.hilbert(target)

    if est.source isa Phase
        s = angle.(hil_s)
    elseif est.source isa Amplitude
        s = abs.(hil_s)
    else
        throw(ArgumentError("est.source must be either Phase or Amplitude instance"))
    end

    if est.target isa Phase
        t = angle.(hil_t)
    elseif est.target isa Amplitude
        t = abs.(hil_t)
    else
        throw(ArgumentError("est.target must be either Phase or Amplitude instance"))
    end

    association(est.est, s, t)
end

function association(est::Hilbert, source, target, cond)
    hil_s = DSP.hilbert(source)
    hil_t = DSP.hilbert(target)
    hil_c = DSP.hilbert(cond)

    if est.source isa Phase
        s = angle.(hil_s)
    elseif est.source isa Amplitude
        s = abs.(hil_s)
    else
        throw(ArgumentError("est.source must be either Phase or Amplitude instance"))
    end

    if est.target isa Phase
        t = angle.(hil_t)
    elseif est.target isa Amplitude
        t = abs.(hil_t)
    else
        throw(ArgumentError("est.target must be either Phase or Amplitude instance"))
    end

    if est.cond isa Phase
        c = angle.(hil_c)
    elseif est.cond isa Amplitude
        c = abs.(hil_c)
    else
        throw(ArgumentError("est.cond must be either Phase or Amplitude instance"))
    end
    
    association(est.est, s, t, c)
end
