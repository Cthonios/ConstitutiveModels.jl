struct UnImplementedMethodError <: Exception
end

abstract type AbstractConstitutive end

# expected interface
function initialize_props end
function initialize_state end
function num_properties end
function num_state_variables end
function property_names end # optional maybe?
function state_variable_names end

# defaults
"""
Default state constructor to just return zeros
"""
function initialize_state(c::AbstractConstitutive)
    return zeros(num_state_variables(c))
end
num_properties(::AbstractConstitutive) = throw(UnImplementedMethodError())
num_state_variables(::AbstractConstitutive) = throw(UnImplementedMethodError())
"""
Return human-readable names for each state variable, in storage order.
Default fallback generates generic names: ["state_1", "state_2", ...].
Models should override this to provide meaningful names.
"""
function state_variable_names(c::AbstractConstitutive)
    return ["state_$i" for i in 1:num_state_variables(c)]
end

"""
for constitutive models which may or may not
be composed of consititutive modules
this interface introduces the expected minimum method inputs of
method(model, props, Z_old, Z_new, Δt)
where model is the model, props is the array of props for this model,
Z_old and Z_new are the old and new state variables and Δt is the time step
"""
abstract type AbstractConstitutiveModel <: AbstractConstitutive end

# minimum interface:
# all models must possess at minimum a density property for now
# this may change in the future to require a density "module"
# for now the first property must be a Lagrangian-frame density value
function density(
    ::AbstractConstitutiveModel,
    props, Z_old, Z_new, Δt, ∇u, θ, args...
)
    return props[1]
end

function dissipation end
function heat_flux end

# below types differentiate between the expected kinematic input for models
# 1. hyperelastic expects the displacement gradient
# 2. hypoelastic expects the velocity gradient and previous stress
# 3. linear elastic expects the linear strain tensor rather than displacement gradient

"""
this adds additional interface expectations
method(model, props, Z_old, Z_new, Δt, ∇u, θ, args...)
"""
abstract type AbstractHyperelasticModel <: AbstractConstitutiveModel end
"""
this add additional interface expections
method(model, props, Z_old, Z_new, Δt, ∇v, θ, σ_old, args...)
"""
abstract type AbstractHypoelasticModel <: AbstractConstitutiveModel end
"""
this adds additional interface expectations
method(model, props, Z_old, Z_new, Δt, ε, θ, args...)
"""
abstract type AbstractLinearElasticModel <: AbstractHyperelasticModel end

# some AD defaults
function entropy(
    model::AbstractHyperelasticModel,
    props, Z_old, Z_new, Δt, ∇u, θ, args...
)
    return -Tensors.gradient(z -> helmholtz_free_energy(
        model, props, Z_old, Z_new, Δt, ∇u, z, args...
    ), θ)
end

function heat_capacity(
    model::AbstractHyperelasticModel,
    props, Z_old, Z_new, Δt, ∇u, θ, args...
)
    return θ * Tensors.gradient(z -> entropy(
        model, props, Z_old, Z_new, Δt, ∇u, z, args...
    ), θ)
end

function material_tangent(
    model::AbstractHyperelasticModel,
    props, Z_old, Z_new, Δt, ∇u, θ, args...
)
    return Tensors.gradient(z -> pk1_stress(
        model, props, Z_old, Z_new, Δt, z, θ, args...
    ), ∇u)
end

function pk1_stress(
    model::AbstractHyperelasticModel,
    props, Z_old, Z_new, Δt, ∇u, θ, args...
)
    return Tensors.gradient(z -> helmholtz_free_energy(
        model, props, Z_old, Z_new, Δt, z, θ, args...
    ), ∇u)
end

function pk1_stress_temperature_modulus(
    model::AbstractHyperelasticModel,
    props, Z_old, Z_new, Δt, ∇u, θ, args...
)
    return Tensors.gradient(z -> pk1_stress(
        model, props, Z_old, Z_new, Δt, ∇u, z, args...
    ), θ)
end

# some default methods, may not always make sense?
function cauchy_stress(
    model::AbstractHyperelasticModel,
    props, Z_old, Z_new, Δt, ∇u, θ, args...
)
    F = ∇u + one(∇u)
    J = det(F)
    P = pk1_stress(model, props, Z_old, Z_new, Δt, ∇u, θ, args...)
    return (1 / J) * dot(P, transpose(F))
end

# linear model defaults
function cauchy_stress(
    model::AbstractLinearElasticModel,
    props, Z_old, Z_new, Δt, ε, θ, args...
)
    return Tensors.gradient(z -> helmholtz_free_energy(
        model, props, Z_old, Z_new, Δt, z, θ, args...
    ), ε)
end

function cauchy_stress_temperature_modulus(
    model::AbstractLinearElasticModel,
    props, Z_old, Z_new, Δt, ε, θ, args...
)
    return Tensors.gradient(z -> cauchy_stress(
        model, props, Z_old, Z_new, Δt, ε, z, args...
    ), θ)
end

function spatial_tangent(
    model::AbstractLinearElasticModel,
    props, Z_old, Z_new, Δt, ε, θ, args...
)
    return Tensors.gradient(z -> cauchy_stress(
        model, props, Z_old, Z_new, Δt, z, θ, args...
    ), ε)
end

# for constitutive "modules"
abstract type AbstractConstitutiveModule <: AbstractConstitutive end

@inline function module_props(
    mod::AbstractConstitutiveModule,
    props,
    start_index::Int
)
    NP = num_properties(mod)
    indices = start_index:(start_index + NP - 1)
    return SVector{NP, eltype(props)}(@views props[indices])
end

@inline function module_props(
    mod::AbstractMaterialSymmetry,
    props,
    start_index::Int
)
    NP = num_properties(mod)
    indices = start_index:(start_index + NP - 1)
    return SVector{NP, eltype(props)}(@views props[indices])
end

# do we really need the below?
abstract type AbstractKinematics end
struct DisplacementGradient <: AbstractKinematics
end
struct LinearStrain <: AbstractKinematics
end
struct VelocityGradient <: AbstractKinematics
end

kinematics(::AbstractHyperelasticModel) = DisplacementGradient()
kinematics(::AbstractHypoelasticModel) = VelocityGradient()
kinematics(::AbstractLinearElasticModel) = LinearStrain()

# some helpers for whether we need certain inputs
requires_temperature_gradient(_) = false
requires_temperature_gradient(::typeof(dissipation)) = true
requires_temperature_gradient(::typeof(heat_flux)) = true