# seperation of phenemenology ownership
# - yield surface owns parameters specific
#   to effective stress calculation
# - isotropic hardening function owns
#   yield stress and other isotropic hardening parameters

include("IsotropicHardening.jl")
include("YieldSurfaces.jl")