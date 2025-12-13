# TODO eventually type correctly
# function linear_strain(∇u::Tensor{2, 3, T, 9}) where T <: Number
function linear_strain(∇u)
    return symmetric(∇u)
end
