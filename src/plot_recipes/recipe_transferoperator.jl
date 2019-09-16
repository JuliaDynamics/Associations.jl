import PerronFrobenius.TransferOperatorRectangularBinning

@recipe function f(r::TransferOperatorRectangularBinning)
    seriestype  :=  :heatmap
    xlabel --> "state #i"
    ylabel --> "state #j"
    r.transfermatrix
end
