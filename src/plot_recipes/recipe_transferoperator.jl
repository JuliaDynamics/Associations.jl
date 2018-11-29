@recipe function f(r::PerronFrobenius.RectangularBinningTransferOperator)
    seriestype  :=  :heatmap
    xlabel --> "state #i"
    ylabel --> "state #j"
    r.transfermatrix
end
