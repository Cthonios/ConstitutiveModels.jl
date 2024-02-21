struct NewtonSolverSettings
  max_iterations::Int32
  absolute_tolerance::Float64
  relative_tolerance::Float64
end

function NewtonSolverSettings()
  return NewtonSolverSettings(Int32(20), 1e-12, 1e-12)
end

function solve_motion(f, x0)
  settings = NewtonSolverSettings()
  R0 = 1.0e6
  x = x0
  for n in 1:settings.max_iterations
    R = f(x)

    if n == 1
      R0 = norm(R)
    end

    if norm(R) < settings.absolute_tolerance
      break
    end

    if norm(R) / R0 < settings.relative_tolerance
      break
    end

    K = ForwardDiff.jacobian(f, x)
    Δx = -K \ R
    x = x + Δx

    if norm(Δx) < settings.absolute_tolerance
      break
    end

    if n == settings.max_iterations
      @warn "Reached maximum Newton iterations. Be careful"
    end
  end

  return x
end

function solve_hardening(f, g, x0)
  settings = NewtonSolverSettings()
  R0 = 1.0e6
  x = x0
  for n in 1:settings.max_iterations
    R = f(x)
    
    if n == 1
      R0 = norm(R)
    end

    if norm(R) < settings.absolute_tolerance
      break
    end

    if norm(R) / R0 < settings.relative_tolerance
      break
    end

    K = g(x)
    Δx = -R / K

    x = x + Δx

    if norm(Δx) < settings.absolute_tolerance
      break
    end

    if n == settings.max_iterations
      @warn "Reached maximum Newton iterations in solve hardening. Be careful"
    end
  end

  return x
end

# function solve_motion(f, x0, state)
#   settings = NewtonSolverSettings()
#   R0 = 1.0e6
#   x = x0
#   for n in 1:settings.max_iterations
#     R = f(x, state)

#     if n == 1
#       R0 = norm(R)
#     end

#     if norm(R) < settings.absolute_tolerance
#       break
#     end

#     if norm(R) / R0 < settings.relative_tolerance
#       break
#     end

#     # K = ForwardDiff.jacobian(f, x)
#     K = ForwardDiff.jacobian(z -> f(z, state), x)
#     Δx = -K \ R
#     x = x + Δx

#     if norm(Δx) < settings.absolute_tolerance
#       break
#     end

#     if n == settings.max_iterations
#       @warn "Reached maximum Newton iterations. Be careful"
#     end
#   end

#   return x
# end
