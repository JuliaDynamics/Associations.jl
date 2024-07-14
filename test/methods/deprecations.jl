@test_throws ArgumentError ExpandingSegment(; libsizes = 10:10:50)
@test_throws ArgumentError RandomVectors(; libsizes = 10:10:50)
@test_throws ArgumentError RandomSegment(; libsizes = 10:10:50)
est =  RandomSegment(CCM(); libsizes = 10:10:50)
x, y = rand(100), rand(100)
@test_throws ArgumentError crossmap(CCM(), est, x, y)

x, y = rand(100), rand(100)
τ = -2
@test_throws ArgumentError crossmap(x, y, 2, τ) isa Real
@test_throws ArgumentError pai(x, y, 2, τ) isa Real
