# From Program 10
using LinearAlgebra

function meshgrid(x, y)
    x = reshape(x, 1, length(x))
    y = reshape(y, length(y), 1)
    X = repeat(x, length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end