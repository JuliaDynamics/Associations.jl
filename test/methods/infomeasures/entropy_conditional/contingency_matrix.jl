# Analytical test from Cover & Thomas textbook.
freqs_yx = [1//8 1//16 1//32 1//32;
    1//16 1//8  1//32 1//32;
    1//16 1//16 1//16 1//16;
    1//4  0//1  0//1  0//1];

freqs_xy = transpose(freqs_yx);
probs_xy = freqs_xy ./ sum(freqs_xy)
c_xy = ContingencyMatrix(probs_xy, freqs_xy)
ce_x_given_y = entropy_conditional(CEShannon(), c_xy) |> Rational
@test ce_x_given_y == 11//8

probs_yx = freqs_yx ./ sum(freqs_yx);
c_yx = ContingencyMatrix(probs_yx, freqs_yx);
ce_y_given_x = entropy_conditional(CEShannon(), c_yx) |> Rational
@test ce_y_given_x == 13//8

# We don't have analytical tests for the other conditional entropies, so just test
# that they successfully compute something.
x = rand(["a", "b", "c"], 100)
y = rand(["hello", "yoyo"], 100)
c = contingency_matrix(x, y)

@test entropy_conditional(CEShannon(), c) >= 0.0
@test entropy_conditional(CETsallisAbe(), c) isa Real
@test entropy_conditional(CETsallisFuruichi(), c) isa Real
