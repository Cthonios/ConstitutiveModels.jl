f(x) = x^3 - 4.
g(x) = sin(x) + x

x0 = 1e-5
x_actual = cbrt(4.)

bracket = (eps(x0), 100.)

solver = ConstitutiveModels.ScalarSolver()

x = ConstitutiveModels.solve(solver, f, x0, bracket)

@test x ≈ x_actual

# test to make sure it gives a Nan
bracket = (2., 100.)
x = ConstitutiveModels.solve(solver, f, x0, bracket)
# @show x
@test isnan(x)

# test a harder function
x0 = 19.
bracket = (-3., 20.)
x = ConstitutiveModels.solve(solver, g, x0, bracket)
@test x ≈ 0.

# test it's differentiable
x0 = 3.
bracket = (eps(x0), 100.)
h(x, a) = x^3 - a
# df = ForwardDiff.derivative(z -> ConstitutiveModels.solve(solver, f, z, bracket), x0)
# @show df

function cube_root(a)
  f_temp(x) = h(x, a)
  bracket = (eps(x0), a)
  ConstitutiveModels.solve(solver, f_temp, 8., bracket)
end

@test cube_root(4.) ≈ cbrt(4.)
@test ForwardDiff.derivative(cube_root, 3.) ≈ 3.0^(-2. / 3.) / 3.

# test with force bisection step

f(x, a) = x^2 - a

function square_root(a)
  f_temp(x) = f(x, a)
  bracket = (eps(x0), a)
  ConstitutiveModels.solve(solver, f_temp, 8., bracket)
end

@test square_root(9.) ≈ 3.

# test when left bracket is solution
f_l(x) = x * (x^2 - 10.)
bracket = (0., 1.)
x0 = 3.
x = ConstitutiveModels.solve(solver, f_l, x0, bracket)
@test x ≈ 0.0

# test when right bracket is solution
f_r(x) = x * (x^2 - 10.)
bracket = (-1., 0.)
x0 = 3.
x = ConstitutiveModels.solve(solver, f_r, x0, bracket)
@test x ≈ 0.0