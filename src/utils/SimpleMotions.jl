# TODO move somewhere else
function _newton(
    f,
    x0;
    atol = 1e-12,
    rtol = 1e-10,
    maxiter = 25,
)
    x = copy(x0)

    for iter in 1:maxiter
        r = f(x)

        if norm(r) ≤ atol + rtol * norm(x)
            return x
        end

        J = ForwardDiff.jacobian(f, x)

        Δx = -(J \ r)

        x .+= Δx

        if norm(Δx) ≤ atol + rtol * norm(x)
            return x
        end
    end

    error("Newton solver failed to converge.")
end

# TODO move somewhere else maybe...
struct MaterialHistory{T, K, M, S, Temp, TG}
    time::T
    kinematics::K
    material_output::M
    state::S
    temperature::Temp
    temperature_gradient::TG
end

abstract type AbstractMotion{
    DGT <: Number, # displacement gradient component type
    VGT <: Number  # velocity gradient component type
} end

function _kinematics(model, motion::AbstractMotion, t, args...)
    if kinematics(model) == DisplacementGradient()
        return displacement_gradient(motion, t, args...)
    elseif kinematics(model) == VelocityGradient()
        return velocity_gradient(motion, t, args...)
    elseif kinematics(model) == LinearStrain()
        return linear_strain(motion, t, args...)
    else
        @assert false "Shouldn't happen"
    end
end

# need to define for each motion
function displacement_gradient end
function velocity_gradient end

"""
$(TYPEDSIGNATURES)
"""
function linear_strain(motion::AbstractMotion, t, args...)
    ∇u = displacement_gradient(motion, t, args...)
    ε = linear_strain(∇u)
    return ε
end

function simulate_material_point(
    mat_func::Function,
    model::AbstractHyperelasticModel,
    props::Dict{String},
    motion::AbstractMotion,
    time_end::Number,
    temp_func = t -> 0.0,
    temp_grad_func = t -> zero(Vec{3, Float64});
    num_steps::Int = 101
)
    # TODO allow for time history that is non-uniform
    Δt = time_end / num_steps

    # initialize props
    props = initialize_props(model, props)

    # initialize state
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)

    times = LinRange(0.0, time_end, num_steps)

    # initialize history outputs
    θ = temp_func(0.0)
    ∇θ = temp_grad_func(0.0)
    if requires_temperature_gradient(mat_func)
        kin_vars = (model, props, Z_old, Z_new, Δt, θ, ∇θ)
        kin = _kinematics(model, motion, 0.0, kin_vars...)
        out = mat_func(model, props, Z_old, Z_new, Δt, kin, θ, ∇θ)
    else
        kin_vars = (model, props, Z_old, Z_new, Δt, 0.0)
        kin = _kinematics(model, motion, 0.0, kin_vars...)
        out = mat_func(model, props, Z_old, Z_new, Δt, kin, θ)
    end

    history_type = MaterialHistory{Float64, typeof(kin), typeof(out), typeof(Z_old), Float64, Vec{3, Float64}}
    history_out = Vector{history_type}(undef, length(times))

    # do loop
    for (n, time) in enumerate(times)
        θ = temp_func(time)
        ∇θ = temp_grad_func(time)
        if requires_temperature_gradient(mat_func)
            kin_vars = (model, props, Z_old, Z_new, Δt, θ, ∇θ)
        else
            kin_vars = (model, props, Z_old, Z_new, Δt, θ)
        end
        kin = _kinematics(model, motion, time, kin_vars...)
        out = mat_func(model, props, Z_old, Z_new, Δt, kin, θ)
        history_out[n] = MaterialHistory(time, kin, out, copy(Z_new), θ, ∇θ)
        Z_old .= Z_new
    end

    return history_out
end

function simulate_material_point(
    mat_func::Function,
    model::AbstractHypoelasticModel,
    props::Dict{String},
    motion::AbstractMotion,
    time_end::Number,
    temp_func = t -> 0.0,
    temp_grad_func = t -> zero(Vec{3, Float64});
    num_steps::Int = 101
)
    # TODO allow for time history that is non-uniform
    Δt = time_end / num_steps

    # initialize props
    props = initialize_props(model, props)

    # initialize state
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)

    times = LinRange(0.0, time_end, num_steps)

    # initialize history outputs
    σ_old = zero(SymmetricTensor{2, 3, Float64, 6})
    θ = temp_func(0.0)
    ∇θ = temp_grad_func(0.0)
    if requires_temperature_gradient(mat_func)
        kin_vars = (model, props, Z_old, Z_new, Δt, σ_old, θ, ∇θ)
        kin = _kinematics(model, motion, 0.0, kin_vars...)
        out = mat_func(model, props, Z_old, Z_new, Δt, σ_old, kin, θ, ∇θ)
    else
        kin_vars = (model, props, Z_old, Z_new, Δt, σ_old, θ)
        kin = _kinematics(model, motion, 0.0, kin_vars...)
        out = mat_func(model, props, Z_old, Z_new, Δt, σ_old, kin, θ)
    end

    history_type = MaterialHistory{Float64, typeof(kin), typeof(out), typeof(Z_old)}
    history_out = Vector{history_type}(undef, length(times))

    # do loop
    for (n, time) in enumerate(times)
        θ = temp_func(time)
        ∇θ = temp_grad_func(time)
        if requires_temperature_gradient(mat_func)
            kin_vars = (model, props, Z_old, Z_new, Δt, σ_old, θ, ∇θ)
            kin = _kinematics(model, motion, time, kin_vars...)
            out = mat_func(model, props, Z_old, Z_new, Δt, σ_old, kin, θ, ∇θ)
        else
            kin_vars = (model, props, Z_old, Z_new, Δt, σ_old, θ)
            kin = _kinematics(model, motion, time, kin_vars...)
            out = mat_func(model, props, Z_old, Z_new, Δt, σ_old, kin, θ)
        end
        history_out[n] = MaterialHistory(time, kin, out, Z_new, θ, ∇θ)
        # update stress
        σ_old = cauchy_stress(model, props, Z_old, Z_new, Δt, σ_old, kin, θ)
        Z_old .= Z_new
    end

    return history_out
end

"""
$(TYPEDEF)
"""
abstract type AbstractConstrainedMotion{DGT, VGT} <: AbstractMotion{DGT, VGT} end

"""
$(TYPEDEF)
"""
abstract type AbstractSimpleMotion{DGT, VGT} <: AbstractMotion{DGT, VGT} end

"""
$(TYPEDSIGNATURES)
"""
function displacement_gradient(motion::AbstractSimpleMotion, t, args...)
    F = deformation_gradient(motion, t, args...)
    I = one(Tensor{2, 3, typeof(t), 9})
    return F - I
end

####################################################
# biaxial strain
####################################################
@doc raw"""
Provides an analytic motion for biaxial strain

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda_1(t) & 0            & 0 \\
0            & \lambda_2(t) & 0 \\
0            & 0            & 1
\end{bmatrix}
```
"""
struct BiaxialStrain{DGT, VGT, F1 <: Function, F2 <: Function} <: AbstractSimpleMotion{DGT, VGT}
    λ1::F1
    λ2::F2

    function BiaxialStrain(λ1::F1, λ2::F2) where {F1 <: Function, F2 <: Function}
        new{Float64, Float64, F1, F2}(λ1, λ2)
    end

    function BiaxialStrain(λ1_dot::T, λ2_dot::T) where T <: Number
        λ1 = t -> λ1_dot * t
        λ2 = t -> λ2_dot * t
        return BiaxialStrain(λ1, λ2)
    end
end

function deformation_gradient(m::BiaxialStrain{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    return Tensor{2, 3, DGT, 9}((
        m.λ1(t), 0,       0,
        0,       m.λ2(t), 0,
        0,       0,       1
    ))
end

function velocity_gradient(m::BiaxialStrain{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    λ1, λ2 = m.λ1(t), m.λ2(t)
    λ1_dot = ForwardDiff.derivative(m.λ1, t)
    λ2_dot = ForwardDiff.derivative(m.λ2, t)
    return Tensor{2, 3, DGT, 9}((
        λ1_dot / λ1, 0,           0,
        0,           λ2_dot / λ2, 0,
        0,           0,           0
    ))
end

####################################################
# isochoric uniaxial stress
####################################################
@doc raw"""
Provides an analytic motion for isochoric uniaxial stress

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda(t) & 0                    & 0 \\
0          & \frac{1}{\lambda}(t) & 0 \\
0          & 0                    & \frac{1}{\lambda}(t)
\end{bmatrix}
```
"""
struct IsochoricUniaxialStress{DGT, VGT, F <: Function} <: AbstractSimpleMotion{DGT, VGT}
    λ::F

    function IsochoricUniaxialStress(λ::F) where F <: Function
        new{Float64, Float64, F}(λ)
    end

    function IsochoricUniaxialStress(λ_dot::T) where T <: Number
        λ = t -> λ_dot * t
        return IsochoricUniaxialStress(λ)
    end
end

function deformation_gradient(m::IsochoricUniaxialStress{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    λ = m.λ(t)
    return Tensor{2, 3, DGT, 9}((
        λ,  0.,           0.,
        0., 1. / sqrt(λ), 0.,
        0., 0.,           1. / sqrt(λ)
    ))
end

function velocity_gradient(m::IsochoricUniaxialStress{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    λ = m.λ(t)
    λ_dot = ForwardDiff.derivative(m.λ, t)
    a = λ_dot / λ
    return Tensor{2, 3, VGT, 9}((
         a,    0,    0,
         0, -a/2,    0,
         0,    0, -a/2
    ))
end

####################################################
# pure shear strain
####################################################
@doc raw"""
Provides an analytic motion for pure shear strain

This is
```math
\mathbf{F} = \frac{1}{2}\begin{bmatrix}
\left(\lambda + \lambda^{-1}\right) & \left(\lambda - \lambda^{-1}\right) & 0 \\
\left(\lambda - \lambda^{-1}\right) & \left(\lambda + \lambda^{-1}\right) & 0 \\
0                                   & 0                                   & 2
\end{bmatrix}
```
"""
struct PureShearStrain{DGT, VGT, F <: Function} <: AbstractSimpleMotion{DGT, VGT}
    λ::F

    function PureShearStrain(λ::F) where F <: Function
        new{Float64, Float64, F}(λ)
    end

    function PureShearStrain(λ_dot::T) where T <: Number
        λ = t -> λ_dot * t
        return PureShearStrain(λ)
    end
end

function deformation_gradient(m::PureShearStrain{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    λ = m.λ(t)
    return 0.5 * Tensor{2, 3, DGT, 9}((
        λ + 1 / λ, λ - 1 / λ, 0,
        λ - 1 / λ, λ + 1 / λ, 0,
        0,         0,         2
    ))
end

function velocity_gradient(m::PureShearStrain{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    λ     = m.λ(t)
    λ_dot = ForwardDiff.derivative(m.λ, t)
    return Tensor{2, 3, VGT, 9}((
        0,         λ_dot / (2λ^2), 0,
        λ_dot / 2, 0,              0,
        0,         0,              0
    ))
end

####################################################
# simple shear
####################################################
@doc raw"""
Provides an analytic motion for simple shear.

This is
```math
\mathbf{F} = \begin{bmatrix}
1 & \gamma & 0 \\
0 & 1      & 0 \\
0 & 0      & 1
\end{bmatrix}
```
"""
struct SimpleShear{DGT, VGT, F <: Function} <: AbstractSimpleMotion{DGT, VGT}
    γ::F

    function SimpleShear(γ::F) where F <: Function
        new{Float64, Float64, F}(γ)
    end

    function SimpleShear(γ_dot::T) where T <: Number
        γ = t -> γ_dot * t
        return SimpleShear(γ)
    end
end

function deformation_gradient(m::SimpleShear{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    return Tensor{2, 3, DGT, 9}((
        1,      0, 0,
        m.γ(t), 1, 0,
        0,      0, 1
    ))
end

function velocity_gradient(m::SimpleShear{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    γ_dot = ForwardDiff.derivative(m.γ, t)
    return Tensor{2, 3, VGT, 9}((
        0,     0, 0,
        γ_dot, 0, 0,
        0,     0, 0
    ))
end

####################################################
# uniaxial strain
####################################################
@doc raw"""
Provides an analytic motion for uniaxial strain

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda(t) & 0 & 0 \\
0          & 1 & 0 \\
0          & 0 & 1
\end{bmatrix}
```

user needs to provide function \lambda(t)
"""
struct UniaxialStrain{DGT, VGT, F <: Function} <: AbstractSimpleMotion{DGT, VGT}
    λ::F

    function UniaxialStrain(λ::F) where F <: Function
        new{Float64, Float64, F}(λ)
    end

    function UniaxialStrain(λ_dot::T) where T <: Number
        λ = t -> λ_dot * t
        return UniaxialStrain(λ)
    end
end

function deformation_gradient(m::UniaxialStrain{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    return Tensor{2, 3, DGT, 9}((
        m.λ(t), 0, 0,
        0,      1, 0,
        0,      0, 1
    ))
end

function velocity_gradient(m::UniaxialStrain{DGT, VGT}, t::T, args...) where {DGT, VGT, T <: Number}
    λ = m.λ(t)
    λ_dot = ForwardDiff.derivative(m.λ, t)
    return Tensor{2, 3, VGT, 9}((
        λ_dot / λ, 0, 0,
        0,         0, 0,
        0,         0, 0
    ))
end

####################################################
# uniaxial stress displacement control
####################################################
struct UniaxialStressDisplacementControl{DGT, VGT, F <: Function} <: AbstractConstrainedMotion{DGT, VGT}
    λ::F

    function UniaxialStressDisplacementControl(λ::F) where F <: Function
        new{Float64, Float64, F}(λ)
    end

    function UniaxialStressDisplacementControl(λ_dot::T) where T <: Number
        λ = t -> λ_dot * t
        return UniaxialStressDisplacementControl(λ)
    end
end

function deformation_gradient(
    motion::UniaxialStressDisplacementControl,
    t,
    model,
    model_inputs...
)
    ∇u = displacement_gradient(motion, t, model, model_inputs...)
    F = ∇u + one(∇u)
    return F
end

function displacement_gradient(
    m::UniaxialStressDisplacementControl{DGT, VGT},
    t::T,
    model, props, Z_old, Z_new, Δt, θ
) where {DGT, VGT, T <: Number}
    function objective(x, t, model, props, Z_old, Z_new, Δt, θ)
        ∇u = Tensor{2, 3, eltype(x), 9}((
            m.λ(t) - 1, 0,    0, 
            0,          x[1], 0, 
            0,          0,    x[1]
        ))
        σ = cauchy_stress(model, props, Z_old, Z_new, Δt, ∇u, θ)
        return [σ[2, 2], σ[3, 3]]
    end
    f = x -> objective(x, t, model, props, Z_old, Z_new, Δt, θ)

    x0 = zeros(1)
    x = _newton(f, x0)[1]
    return Tensor{2, 3, DGT, 9}((
        m.λ(t) - 1, 0,    0, 
        0,          x[1], 0, 
        0,          0,    x[1]
    ))
end

function linear_strain(
    motion::UniaxialStressDisplacementControl{DGT, VGT},
    t::T,
    model, props, Z_old, Z_new, Δt, θ
) where {DGT, VGT, T <: Number}

    # if λ ≈ 1
    #     return zero(SymmetricTensor{2, 3, T, 6})
    # end

    function objective(x, t, model, props, Z_old, Z_new, Δt, θ)
        ε = SymmetricTensor{2, 3, eltype(x), 6}((
            motion.λ(t) - 1, 0, 0, x[1], 0, x[1]
        ))
        σ = cauchy_stress(model, props, Z_old, Z_new, Δt, ε, θ)
        return [σ[2, 2], σ[3, 3]]
    end
    f = x -> objective(x, t, model, props, Z_old, Z_new, Δt, θ)

    x0 = zeros(1)
    x = _newton(f, x0)[1]
    return SymmetricTensor{2, 3, DGT, 6}((
        motion.λ(t) - 1, 0, 0, x, 0, x
    ))
end
