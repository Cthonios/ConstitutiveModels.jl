using ConstitutiveModels
using LinearAlgebra
using Plots
using StaticArrays
using Tensors

inputs = Dict(
  "Young's modulus"           => 70e9,
  "Poisson's ratio"           => 0.3,
  "yield surface"             => "VonMisesYieldSurface",
  "yield stress"              => 200.0e6,
  "isotropic hardening model" => "LinearIsotropicHardening",
  "hardening modulus"         => 200.0e6
)

# model = TrescaYieldSurface(
model = VonMisesYieldSurface(
    LinearIsotropicHardening()
)
props = initialize_props(model, inputs)

# Set up Lode angle θ from -π/6 to π/6 (the full π-plane)
function deviatoric_direction(θ::Float64)
    # Construct 2 orthonormal deviatoric directions
    e1 = SVector(2/√6, -1/√6, -1/√6, 0.0, 0.0, 0.0)  # axial
    e2 = SVector(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)       # shear in 1–2 plane
    en = normalize(cos(θ) * e1 + sin(θ) * e2)
    return SymmetricTensor{2, 3, Float64, 6}((
        en[1], en[4], en[6], en[2], en[3], en[5]
    ))
end

# # Von Mises yield function: √(3/2 s_ij s_ij)
function von_mises(s)
    # return sqrt(3/2 * dot(s, s))  # assumes s is already deviatoric
    return sqrt(3/2) * norm(dev(s))
end

function tresca(s)
    vals = eigvals(s)
    return max(
        abs(vals[1] - vals[2]),
        abs(vals[2] - vals[3]),
        abs(vals[3] - vals[1])
    )
end

# # Compute yield radius r(θ) such that f(r * ŝ) = σ₀ (usually 1.0)
# function yield_radius(eqps; n=3000, σ₀=1.0)
#     # θ = range(-π/6, π/6; length=n)
#     θ = range(0, 2π, length=n)
#     r_vals = Float64[]

#     for θi in θ
#         s = deviatoric_direction(θi)
#         # f_val = yield_function(s)
#         # f_val = ConstitutiveModels.radius(model.isotropic_hardening, props, eqps)
#         f_val = ConstitutiveModels.effective_stress(model, s)^2
#         r = σ₀ / f_val
#         push!(r_vals, r)
#     end

#     x = r_vals .* cos.(θ)
#     y = r_vals .* sin.(θ)
#     return x, y, θ
# end

# # Plotting von Mises yield surface
# x_vm, y_vm, _ = yield_radius(0.0)

# p = plot(x_vm, y_vm)#, label="Von Mises", linewidth=2, aspect_ratio=1, legend=:bottomleft)
# xlabel!("π-plane x")
# ylabel!("π-plane y")
# title!("Yield Surface in π-plane")

# eqps_vals = range(0., 0.1, length=5)
# for eqps in eqps_vals
#     x_vm, y_vm, _ = yield_radius(eqps)
#     plot!(p, x_vm, y_vm)#, label="Von Mises", linewidth=2, aspect_ratio=1, legend=:bottomleft)
# end
# p

# Convert principal stresses to deviatoric form
function deviatoric_stress(σ::SVector{3})
    p = sum(σ) / 3
    return σ .- p
end

# Project deviatoric stress in 3D to 2D π-plane coordinates
function pi_plane_coords(σ_dev::SVector{3})
    # Use standard π-plane orthonormal basis:
    # e_ξ = (1, -1, 0)/√2
    # e_η = (1, 1, -2)/√6
    ξ = dot(σ_dev, SVector(1.0, -1.0, 0.0)) / √2
    η = dot(σ_dev, SVector(1.0, 1.0, -2.0)) / √6
    return (-ξ, -η)
end

# Principal directions as points
σ1 = SVector(1.0, 0.0, 0.0)
σ2 = SVector(0.0, 1.0, 0.0)
σ3 = SVector(0.0, 0.0, 1.0)

# Deviatoric parts
s1 = deviatoric_stress(σ1)
s2 = deviatoric_stress(σ2)
s3 = deviatoric_stress(σ3)

# Project to 2D
x1, y1 = pi_plane_coords(s1)
x2, y2 = pi_plane_coords(s2)
x3, y3 = pi_plane_coords(s3)

# Plot arrows from origin
p = plot(arrow=(:closed, 0.2), aspect_ratio=1, legend=false, size=(500, 500))
quiver!([0.0], [0.0], quiver=([x1], [y1]), label="σ₁", color=:red)
quiver!([0.0], [0.0], quiver=([x2], [y2]), label="σ₂", color=:blue)
quiver!([0.0], [0.0], quiver=([x3], [y3]), label="σ₃", color=:green)

# Decorations
xlabel!("π-plane ξ")
ylabel!("π-plane η")
title!("Principal Stress Axes Projected onto π-plane")

θ = range(0, 2π, length=30000)
r_vals_vm = Float64[]
r_vals_tr = Float64[]

for θi in θ
    s = deviatoric_direction(θi)
    f_val_vm = von_mises(s)
    f_val_tr = tresca(s)
    r_vm = 1. / f_val_vm
    r_tr = 1. / f_val_tr
    push!(r_vals_vm, r_vm)
    push!(r_vals_tr, r_tr)
end
x = r_vals_vm .* cos.(θ)
y = r_vals_vm .* sin.(θ)
plot!(p, x, y)
x = r_vals_tr .* cos.(θ)
y = r_vals_tr .* sin.(θ)
plot!(p, x, y)