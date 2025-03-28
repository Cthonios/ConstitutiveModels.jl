struct ScalarSolver
  max_iters::Int
  x_tol::Float64
  r_tol::Float64
end

function ScalarSolver()
  return ScalarSolver(50, 1e-13, 0)
end

function clip(x, l, h)
  if x < l
    return l
  elseif x > h
    return h
  else
    return x
  end
end

function solve(solver::ScalarSolver, f, x0, bracket)
  converged = false
  df = z -> ForwardDiff.derivative(f, z)
  fl, fh = f(bracket[1]), f(bracket[2])

  x0 = clip(x0, bracket...)

  # check that root is bracketed
  x0 = ifelse(fl * fh < zero(eltype(fl)), x0, NaN)

  # Check if either bracket is a root
  left_bracket_is_solution = (fl == 0.0)
  x0 = ifelse(left_bracket_is_solution, bracket[1], x0)
  converged = ifelse(left_bracket_is_solution, true, converged)

  right_bracked_is_solution = (fh == 0.0)
  x0 = ifelse(right_bracked_is_solution, bracket[2], x0)
  converged = ifelse(right_bracked_is_solution, true, converged)

  # ORIENT THE SEARCH SO THAT F(XL) < 0.
  xl, xh = ifelse(fl < 0, (bracket[1], bracket[2]), (bracket[2], bracket[1]))
  
  # INITIALIZE THE ''STEP SIZE BEFORE LAST'', AND THE LAST STEP
  dx_old = abs(bracket[2] - bracket[1])
  dx = dx_old

  F, DF = f(x0), df(x0)

  # f_calls = 3

  x = x0

  i = 1
  while !converged && i < solver.max_iters
    newton_out_of_range = ((x - xh) * DF - F) * ((x - xl) * DF - F) > 0
    newton_decrease_slowly = abs(2 * F) > abs(dx_old * DF)
    dx_old = dx

    if newton_out_of_range || newton_decrease_slowly
      # take bisection step
      dx = 0.5 * (xh - xl)
      x = xl + dx
      converged = (x == xl)
    else
      # take newton step
      dx = -F / DF
      temp = x
      x = x + dx
      converged = (x == temp)
    end

    F, DF = f(x), df(x)

    # maintaint he bracket on the root
    xl, xh = ifelse(F < 0, (x, xh), (xl, x))
    i = i + 1
    coverged = converged || abs(dx) < solver.x_tol || abs(F) < solver.r_tol
  end

  x = ifelse(converged, x, NaN)
  return x
end
