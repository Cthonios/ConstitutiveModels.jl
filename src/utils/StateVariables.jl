function state_variable_names(::SymmetricTensor{2, 3, T, 6}, name::String) where T <: Number
    exts = ["_xx", "_xy", "_xz", "_yy", "_yz", "_zz"]
    return map(x -> name * x, exts)
end

function state_variable_names(::Tensor{2, 3, T, 9}, name::String) where T <: Number
    exts = [
        "_xx", "_yx", "_zx",
        "_xy", "_yy", "_zy",
        "_xz", "_yz", "_zz"
    ]
    return map(x -> name * x, exts)
end
