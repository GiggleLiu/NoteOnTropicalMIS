using .CUDA

function onehotmask(A::CuArray{T}, X::CuArray{T}) where T
    CuArray(onehotmask(Array(A), Array(X)))
end