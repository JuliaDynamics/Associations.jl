# Overview

The association measure API consists of several distinct subcomponents that all have 
one goal: to quantify the "association" between input datasets in some way. Precisely
what is meant by "association" depends on the measure, and precisely what is meant 
by "quantify" depends on the *estimator* of that measure. You'll find more 
info under the following subheadings.

- The [information API](@ref information_api).
- The [cross mapping API](@ref cross_mapping_api).
- The [correlation API](@ref correlation_api).

Association measures are useful in their own right, but their utility come into 
play when using them in the context of [independence testing](@ref independence_testing).