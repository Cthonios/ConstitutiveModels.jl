# the below code was copied from github.com/julialang/julia
# from the file base/special/log.jl to disable the complex
# domain throw error and return NaN instead so that
# solvers that can handle NaNs (such as a trust region solver)
# can move forward with out this throw killing the simulations
#
# it also gives app devs the oppurtuntiy to handle errors
# such as this at a higher level to e.g. report things
# like inverted element or something like that.

# @noinline safe_log(x::T) where T <: Number = safe_log(x)
@noinline safe_log(x::T) where T <: Number = _safe_log(x, Val(:ℯ))

@inline function _safe_log(x::Float64, base)
    if x > 0.0
        x == Inf && return x

        # Step 2
        if 0.9394130628134757 < x < 1.0644944589178595
            f = x-1.0
            return Base.Math.log_proc2(f, base)
        end

        # Step 3
        xu = reinterpret(UInt64,x)
        # `% Int` rather than `Int(...)` to preserve `:nothrow` (the shifted value fits in 11 bits)
        m = ((xu >> 52) % Int) & 0x07ff
        if m == 0 # x is subnormal
            x *= 1.8014398509481984e16 # 0x1p54, normalise significand
            xu = reinterpret(UInt64,x)
            m = ((xu >> 52) % Int) & 0x07ff - 54
        end
        m -= 1023
        y = reinterpret(Float64,(xu & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000)

        mf = Float64(m)
        F = (y + 3.5184372088832e13) - 3.5184372088832e13 # 0x1p-7*round(0x1p7*y)
        f = y-F

        return Base.Math.log_proc1(y,mf,F,f,base)
    elseif x == 0.0
        return -Inf
    elseif isnan(x)
        return NaN
    else
        # do not throw a complex domain error
        # instead throw a NaN so we can handle this later
        # throw_complex_domainerror(func, x)
        return NaN
    end
end

# Tensors.jl implementation
# @implement_gradient safe_log inv
@inline function safe_log_gradient(x::T) where {T <: Number}
    return safe_log(x), inv(x)
end

@implement_gradient safe_log safe_log_gradient

