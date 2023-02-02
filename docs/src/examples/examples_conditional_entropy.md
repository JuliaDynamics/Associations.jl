# [Entropy](@id examples_condentropy)

## Discrete: example from Cover & Thomas

This is essentially example 2.2.1 in Cover & Thomas (2006), where they use the following
contingency table as an example. We'll take their example and manually construct
a [`ContingencyMatrix`](@ref) that we can use to compute the conditional entropy.
The [`ContingencyMatrix`](@ref) constructor takes the probabilities as the
first argument and the raw frequencies as the second argument.
Note also that Julia is column-major, so we need to transpose their example. Then their
`X` is in the first dimension of our contingency matrix (along columns) and their `Y` is
our second dimension (rows).

```@example ce_contingency_table
using CausalityTools
freqs_yx = [1//8 1//16 1//32 1//32; 
    1//16 1//8  1//32 1//32;
    1//16 1//16 1//16 1//16; 
    1//4  0//1  0//1  0//1];
freqs_xy = transpose(freqs_yx);
probs_xy = freqs_xy ./ sum(freqs_xy)
c_xy = ContingencyMatrix(probs_xy, freqs_xy)
```

The marginal distribution for `x` (first dimension) is

```@example ce_contingency_table
probabilities(c_xy, 1)
```

The marginal distribution for `y` (second dimension) is

```@example ce_contingency_table
probabilities(c_xy, 2)
```

And the Shannon conditional entropy ``H^S(X | Y)``

```@example ce_contingency_table
ce_x_given_y = entropy_conditional(CEShannon(), c_xy) |> Rational
```

This is the same as in their example. Hooray! To compute ``H^S(Y | X)``, we just need to
flip the contingency matrix.

```@example ce_contingency_table
probs_yx = freqs_yx ./ sum(freqs_yx);
c_yx = ContingencyMatrix(probs_yx, freqs_yx);
ce_y_given_x = entropy_conditional(CEShannon(), c_yx) |> Rational
```
