abstract type AbstractIsotropicHardening <: AbstractConstitutiveModule end
num_state_variables(::AbstractIsotropicHardening) = 1
yield_surface_radius(model::AbstractIsotropicHardening, props, eqps) = sqrt(2. / 3.) * hardening_gradient(model, props, eqps)

# expected inteface below
function hardening_energy end
function hardening_gradient end
function hardening_hessian end

struct NoIsotropicHardening <: AbstractIsotropicHardening
end

function initialize_props(::NoIsotropicHardening, inputs::Dict{String})
    return [get_property(inputs, "yield stress")]
end

num_properties(::NoIsotropicHardening) = 1

hardening_energy(::NoIsotropicHardening, props, eqps) = props[1] * eqps
hardening_gradient(::NoIsotropicHardening, props, eqps) = props[1]
hardening_hessian(::NoIsotropicHardening, props, eqps) = zero(typeof(eqps))

struct LinearIsotropicHardening <: AbstractIsotropicHardening
end

function initialize_props(::LinearIsotropicHardening, inputs::Dict{String})
    return [
        get_property(inputs, "yield stress"),
        get_property(inputs, "hardening modulus")
    ]
end

num_properties(::LinearIsotropicHardening) = 2

hardening_energy(::LinearIsotropicHardening, props, eqps) = props[1] * eqps + 0.5 * props[2] * eqps * eqps
hardening_gradient(::LinearIsotropicHardening, props, eqps) = props[1] + props[2] * eqps
hardening_hessian(::LinearIsotropicHardening, props, eqps) = props[2]

struct PowerLawIsotropicHardening <: AbstractIsotropicHardening
end

function initialize_props(::PowerLawIsotropicHardening, inputs::Dict{String})
    return [
        get_property(inputs, "yield stress"),         # σ_y
        get_property(inputs, "hardening coefficient"),# K
        get_property(inputs, "hardening exponent"),   # n
        get_property(inputs, "luders strain")         # ε_L
    ]
end

num_properties(::PowerLawIsotropicHardening) = 4

function hardening_energy(::PowerLawIsotropicHardening, props, eqps)
    σy = props[1]
    K  = props[2]
    n  = props[3]
    εL = props[4]

    if eqps <= εL
        return σy * eqps
    else
        Δ = eqps - εL
        return σy * eqps +
               K / (n + 1) * Δ^(n + 1)
    end
end

function hardening_gradient(::PowerLawIsotropicHardening, props, eqps)
    σy = props[1]
    K  = props[2]
    n  = props[3]
    εL = props[4]

    if eqps <= εL
        return σy
    else
        Δ = eqps - εL
        return σy + K * Δ^n
    end
end

function hardening_hessian(::PowerLawIsotropicHardening, props, eqps)
    K  = props[2]
    n  = props[3]
    εL = props[4]

    if eqps <= εL
        return zero(eqps)
    else
        Δ = eqps - εL
        return K * n * Δ^(n - 1)
    end
end

struct Voce <: AbstractIsotropicHardening
end

function initialize_props(::Voce, inputs::Dict{String})
    return [
        get_property(inputs, "yield stress"),
        get_property(inputs, "hardening modulus"),
        get_property(inputs, "saturation stress"),
        get_property(inputs, "hardening exponent")
    ]
end

num_properties(::Voce) = 4

# hardening_energy(::Voce, props, eqps) = props[1] * eqps + 0.5 * props[2] * eqps * eqps + props[3] * eqps - 
# hardening_gradient(::Voce, props, eqps) = props[1] + props[2] * eqps
# hardening_hessian(::Voce, props, eqps) = props[2]

function hardening_energy(::Voce, props, eqps)
    σ_y = props[1]
    H   = props[2]
    A   = props[3]
    n   = props[4]

    return σ_y * eqps +
           0.5 * H * eqps^2 +
           A * eqps -
           (A / n) * (1 - exp(-n * eqps))
end

function hardening_gradient(::Voce, props, eqps)
    σ_y = props[1]
    H   = props[2]
    A   = props[3]
    n   = props[4]

    return σ_y +
           H * eqps +
           A * (1 - exp(-n * eqps))
end

function hardening_hessian(::Voce, props, eqps)
    H = props[2]
    A = props[3]
    n = props[4]

    return H + A * n * exp(-n * eqps)
end
