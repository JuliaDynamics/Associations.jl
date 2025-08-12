# Check if the difference is within a certain threshold percentage. Used to check 
# agreement between `ValueBinning` and `TransferOperator` estimation.
# `agreement_threshold` is a percentage indicating the maximal percentage 
# discrepancy between the two values.
function in_agreement(val1, val2; agreement_threshold=0.05)
    largest_magnitude = max(abs(val1), abs(val2))
    if largest_magnitude == 0
        return val1 == val2
    else
        return abs(val1 - val2) / largest_magnitude <= agreement_threshold
    end
end