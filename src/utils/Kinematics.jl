# TODO eventually type correctly
# function linear_strain(∇u::Tensor{2, 3, T, 9}) where T <: Number
function linear_strain(∇u)
    return symmetric(∇u)
end

function _sort_principal_stretches(a::T, b::T, c::T) where T<:Real
    vals = (a, b, c)
    inds = (1, 2, 3)

    # Compare first two
    if vals[1] > vals[2]
        vals = (vals[2], vals[1], vals[3])
        inds = (inds[2], inds[1], inds[3])
    end
    # Compare last two
    if vals[2] > vals[3]
        vals = (vals[1], vals[3], vals[2])
        inds = (inds[1], inds[3], inds[2])
    end
    # Compare first two again
    if vals[1] > vals[2]
        vals = (vals[2], vals[1], vals[3])
        inds = (inds[2], inds[1], inds[3])
    end

    return vals, inds
end

function principal_stretchs(A::SymmetricTensor{2, 3, T, 6}) where T <: Real
    # invariants of A
    I1 = tr(A)
    I2 = 0.5*(I1^2 - tr(A ⋅ A))
    I3 = det(A)

    # cubic: x^3 - I1 x^2 + I2 x - I3 = 0
    b = -I1
    c = I2
    d = -I3

    Δ0 = b^2 - 3*c
    Δ1 = 2*b^3 - 9*b*c + 27*d

    # safe computation of theta
    θ = acos(clamp(Δ1 / (2 * sqrt(Δ0^3)), -1.0, 1.0))

    # three real roots (Cardano trig form)
    x1 = -(b/3 + 2*sqrt(Δ0)/3 * cos(θ/3))
    x2 = -(b/3 + 2*sqrt(Δ0)/3 * cos((θ + 2π)/3))
    x3 = -(b/3 + 2*sqrt(Δ0)/3 * cos((θ + 4π)/3))

    # principal stretches
    (x1, x2, x3), _ = _sort_principal_stretches(x1, x2, x3)
    return Vec{3, T}((sqrt(x1), sqrt(x2), sqrt(x3)))
end
