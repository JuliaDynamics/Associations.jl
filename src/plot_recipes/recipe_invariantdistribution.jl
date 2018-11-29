@recipe function f(id::PerronFrobenius.InvariantDistribution)
    seriestype  :=  :bar
    xlabel --> "state # (i)"
    ylabel --> "invariant density"
    linecolor --> :black
    fillcolor --> :black
    fillalpha --> 0.5
    legend --> :none
    id.dist
end
