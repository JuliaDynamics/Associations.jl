abstract type AbstractDynamicalSystemResampling end

"""
    NoResampling

Indicates that no resampling should be performed.
"""
struct NoResampling <: AbstractDynamicalSystemResampling end 


export AbstractDynamicalSystemResampling, NoResampling