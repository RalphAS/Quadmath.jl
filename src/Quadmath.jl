module Quadmath
using Requires

export Float128, ComplexF128

import Base: (*), +, -, /,  <, <=, ==, ^, convert,
          reinterpret, sign_mask, exponent_mask, exponent_one, exponent_half,
          significand_mask, exponent, significand,
          promote_rule, widen,
          string, print, show, parse,
          acos, acosh, asin, asinh, atan, atanh, cosh, cos,
          exp, expm1, log, log2, log10, log1p, sin, sinh, sqrt,
          tan, tanh,
          ceil, floor, trunc, round, fma,
          copysign, flipsign, max, min, hypot, abs,
          ldexp, frexp, modf, nextfloat, eps,
          isinf, isnan, isfinite, isinteger,
          floatmin, floatmax, precision, signbit,
          Int32, Int64, Float64, BigFloat, BigInt

if Sys.isapple()
    const quadoplib = "libquadmath.0"
    const libquadmath = "libquadmath.0"
elseif Sys.isunix()
    const quadoplib = "libgcc_s.so.1"
    const libquadmath = "libquadmath.so.0"
elseif Sys.iswindows()
    const quadoplib = "libgcc_s_seh-1.dll"
    const libquadmath = "libquadmath-0.dll"
end

@static if Sys.isunix()
    # we use this slightly cumbersome definition to ensure that the value is passed
    # on the xmm registers, matching the x86_64 ABI for __float128.
    const Cfloat128 = NTuple{2,VecElement{Float64}}

    struct Float128 <: AbstractFloat
        data::Cfloat128
    end
    convert(::Type{Float128}, x::Number) = Float128(x)

    const ComplexF128 = Complex{Float128}

    Base.cconvert(::Type{Cfloat128}, x::Float128) = x.data


    # reinterpret
    function reinterpret(::Type{UInt128}, x::Float128)
        hi = reinterpret(UInt64, x.data[2].value)
        lo = reinterpret(UInt64, x.data[1].value)
        UInt128(hi) << 64 | lo
    end
    function reinterpret(::Type{Float128}, x::UInt128)
        fhi = reinterpret(Float64, (x >> 64) % UInt64)
        flo = reinterpret(Float64, x % UInt64)
        Float128((VecElement(flo), VecElement(fhi)))
    end
    reinterpret(::Type{Unsigned}, x::Float128) = reinterpret(UInt128, x)
    reinterpret(::Type{Signed}, x::Float128) = reinterpret(Int128, x)

    reinterpret(::Type{Int128}, x::Float128) =
        reinterpret(Int128, reinterpret(UInt128, x))
    reinterpret(::Type{Float128}, x::Int128) =
        reinterpret(Float128, reinterpret(UInt128, x))

elseif Sys.iswindows()
    primitive type Float128 <: AbstractFloat 128
    end
    const Cfloat128 = Float128
end

function __init__()
    @require SpecialFunctions="276daf66-3868-5448-9aa4-cd146d93841b" begin
        import .SpecialFunctions

        SpecialFunctions.erf(x::Float128) = Float128(ccall((:erfq, libquadmath), Cfloat128, (Cfloat128, ), x))
        SpecialFunctions.erfc(x::Float128) = Float128(ccall((:erfcq, libquadmath), Cfloat128, (Cfloat128, ), x))

        SpecialFunctions.besselj0(x::Float128) = Float128(ccall((:j0q, libquadmath), Cfloat128, (Cfloat128, ), x))
        SpecialFunctions.besselj1(x::Float128) = Float128(ccall((:j1q, libquadmath), Cfloat128, (Cfloat128, ), x))
        SpecialFunctions.bessely0(x::Float128) = Float128(ccall((:y0q, libquadmath), Cfloat128, (Cfloat128, ), x))
        SpecialFunctions.bessely1(x::Float128) = Float128(ccall((:y1q, libquadmath), Cfloat128, (Cfloat128, ), x))
        SpecialFunctions.besselj(n::Cint, x::Float128) = Float128(ccall((:jnq, libquadmath), Cfloat128, (Cint, Cfloat128), n, x))
        SpecialFunctions.bessely(n::Cint, x::Float128) = Float128(ccall((:ynq, libquadmath), Cfloat128, (Cint, Cfloat128), n, x))

        SpecialFunctions.gamma(x::Float128) = Float128(ccall((:tgammaq, libquadmath), Cfloat128, (Cfloat128, ), x))
        SpecialFunctions.lgamma(x::Float128) = Float128(ccall((:lgammaq, libquadmath), Cfloat128, (Cfloat128, ), x))
    end
end


sign_mask(::Type{Float128}) =        0x8000_0000_0000_0000_0000_0000_0000_0000
exponent_mask(::Type{Float128}) =    0x7fff_0000_0000_0000_0000_0000_0000_0000
exponent_one(::Type{Float128}) =     0x3fff_0000_0000_0000_0000_0000_0000_0000
exponent_half(::Type{Float128}) =    0x3ffe_0000_0000_0000_0000_0000_0000_0000
significand_mask(::Type{Float128}) = 0x0000_ffff_ffff_ffff_ffff_ffff_ffff_ffff

fpinttype(::Type{Float128}) = UInt128

# conversion
Float128(x::Float128) = x

## Float64
Float128(x::Float64) =
    Float128(ccall((:__extenddftf2, quadoplib), Cfloat128, (Cdouble,), x))
Float64(x::Float128) =
    ccall((:__trunctfdf2, quadoplib), Cdouble, (Cfloat128,), x)

Int32(x::Float128) =
    ccall((:__fixtfsi, quadoplib), Int32, (Cfloat128,), x)
Float128(x::Int32) =
    Float128(ccall((:__floatsitf, quadoplib), Cfloat128, (Int32,), x))

Float128(x::UInt32) =
    Float128(ccall((:__floatunsitf, quadoplib), Cfloat128, (UInt32,), x))

Int64(x::Float128) =
    ccall((:__fixtfdi, quadoplib), Int64, (Cfloat128,), x)
Float128(x::Int64) =
    Float128(ccall((:__floatditf, quadoplib), Cfloat128, (Int64,), x))

Float128(x::Rational{T}) where T = Float128(numerator(x))/Float128(denominator(x))

# comparison

(==)(x::Float128, y::Float128) =
    ccall((:__eqtf2,quadoplib), Cint, (Cfloat128,Cfloat128), x, y) == 0

(<)(x::Float128, y::Float128) =
    ccall((:__letf2,quadoplib), Cint, (Cfloat128,Cfloat128), x, y) == -1

(<=)(x::Float128, y::Float128) =
    ccall((:__letf2,quadoplib), Cint, (Cfloat128,Cfloat128), x, y) <= 0

# arithmetic

(+)(x::Float128, y::Float128) =
    Float128(ccall((:__addtf3,quadoplib), Cfloat128, (Cfloat128,Cfloat128), x, y))
(-)(x::Float128, y::Float128) =
    Float128(ccall((:__subtf3,quadoplib), Cfloat128, (Cfloat128,Cfloat128), x, y))
(*)(x::Float128, y::Float128) =
    Float128(ccall((:__multf3,quadoplib), Cfloat128, (Cfloat128,Cfloat128), x, y))
(/)(x::Float128, y::Float128) =
    Float128(ccall((:__divtf3,quadoplib), Cfloat128, (Cfloat128,Cfloat128), x, y))
(-)(x::Float128) =
    Float128(ccall((:__negtf2,quadoplib), Cfloat128, (Cfloat128,), x))
(^)(x::Float128, y::Float128) =
    Float128(ccall((:powq, libquadmath), Cfloat128, (Cfloat128,Cfloat128), x, y))
# math

## one argument
for f in (:acos, :acosh, :asin, :asinh, :atan, :atanh, :cosh, :cos,
          :erf, :erfc, :exp, :expm1, :log, :log2, :log10, :log1p,
          :sin, :sinh, :sqrt, :tan, :tanh,
          :ceil, :floor, :trunc, )
    @eval function $f(x::Float128)
        Float128(ccall(($(string(f,:q)), libquadmath), Cfloat128, (Cfloat128, ), x))
    end
end
for (f,fc) in (:abs => :fabs,
               :round => :rint,)
    @eval function $f(x::Float128)
        Float128(ccall(($(string(fc,:q)), libquadmath), Cfloat128, (Cfloat128, ), x))
    end
end

## two argument
for f in (:copysign, :hypot, )
    @eval function $f(x::Float128, y::Float128)
       Float128(ccall(($(string(f,:q)), libquadmath), Cfloat128, (Cfloat128, Cfloat128), x, y))
    end
end

flipsign(x::Float128, y::Float128) = signbit(y) ? -x : x

function atan(x::Float128, y::Float128)
    Float128(ccall((:atan2q, libquadmath), Cfloat128, (Cfloat128, Cfloat128), x, y))
end

## misc
fma(x::Float128, y::Float128, z::Float128) =
    Float128(ccall((:fmaq,libquadmath), Cfloat128, (Cfloat128, Cfloat128, Cfloat128), x, y, z))
    
isnan(x::Float128) =
    0 != ccall((:isnanq,libquadmath), Cint, (Cfloat128, ), x)
isinf(x::Float128) =
    0 != ccall((:isinfq,libquadmath), Cint, (Cfloat128, ), x)
isfinite(x::Float128) =
    0 != ccall((:finiteq,libquadmath), Cint, (Cfloat128, ), x)
isinteger(x::Float128) = isfinite(x) && x === trunc(x)

signbit(x::Float128) = signbit(reinterpret(Int128, x))
precision(::Type{Float128}) = 113

eps(::Type{Float128}) = reinterpret(Float128, 0x3f8f_0000_0000_0000_0000_0000_0000_0000)
floatmin(::Type{Float128}) = reinterpret(Float128, 0x0001_0000_0000_0000_0000_0000_0000_0000)
floatmax(::Type{Float128}) = reinterpret(Float128, 0x7ffe_ffff_ffff_ffff_ffff_ffff_ffff_ffff)

ldexp(x::Float128, n::Cint) =
    Float128(ccall((:ldexpq, libquadmath), Cfloat128, (Cfloat128, Cint), x, n))
ldexp(x::Float128, n::Integer) =
    ldexp(x, clamp(n, typemin(Cint), typemax(Cint)) % Cint)

function frexp(x::Float128)
    r = Ref{Cint}()
    y = Float128(ccall((:frexpq, libquadmath), Cfloat128, (Cfloat128, Ptr{Cint}), x, r))
    return y, Int(r[])
end

function modf(x::Float128)
    isinf(x) && return (zero(Float128), x)
    ipart = trunc(x)
    fpart = x - ipart
    return fpart, ipart
end

significand(x::Float128) = frexp(x)[1] * 2
function exponent(x::Float128)
     !isfinite(x) && throw(DomainError("Cannot be NaN or Inf."))
     abs(x) > 0 && return frexp(x)[2] - 1
     throw(DomainError("Cannot be subnormal converted to 0."))
end

function nextfloat(f::Float128, d::Integer)
    F = typeof(f)
    fumax = reinterpret(Unsigned, F(Inf))
    U = typeof(fumax)

    isnan(f) && return f
    fi = reinterpret(Signed, f)
    fneg = fi < 0
    fu = unsigned(fi & typemax(fi))

    dneg = d < 0
    da = Base.uabs(d)
    if da > typemax(U)
        fneg = dneg
        fu = fumax
    else
        du = da % U
        if fneg ⊻ dneg
            if du > fu
                fu = min(fumax, du - fu)
                fneg = !fneg
            else
                fu = fu - du
            end
        else
            if fumax - fu < du
                fu = fumax
            else
                fu = fu + du
            end
        end
    end
    if fneg
        fu |= sign_mask(F)
    end
    reinterpret(F, fu)
end

Float128(::Irrational{:π}) =  reinterpret(Float128, 0x4000921fb54442d18469898cc51701b8)
Float128(::Irrational{:e}) =  reinterpret(Float128, 0x40005bf0a8b1457695355fb8ac404e7a)

function BigFloat(x::Float128; precision=precision(BigFloat))
    if !isfinite(x) || iszero(x)
        @static if VERSION < v"1.1"
            return BigFloat(Float64(x), precision)
        else
            return BigFloat(Float64(x), precision=precision)
        end
    end

    @static if VERSION < v"1.1"
        b = setprecision(BigFloat, max(precision,113)) do
            BigFloat()
        end
    else
        b = BigFloat(precision=max(precision,113))
    end

    y, k = frexp(x)
    b.exp = Clong(k)
    b.sign = signbit(x) ? Cint(-1) : Cint(1)
    u = (reinterpret(UInt128, y) << 15) | 0x8000_0000_0000_0000_0000_0000_0000_0000
    i = cld(precision, sizeof(Culong)*8)
    while u != 0
        w = (u >> (128-sizeof(Culong)*8)) % Culong
        unsafe_store!(b.d, w, i)
        i -= 1
        u <<= sizeof(Culong)*8
    end
    # set remaining bits to zero
    while i > 0
        unsafe_store!(b.d, zero(Culong), i)
        i -= 1
    end

    if precision < 113
        @static if VERSION < v"1.1"
            b2 = setprecision(BigFloat, precision) do
                BigFloat()
            end
            ccall((:mpfr_set, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
                  b2, b, Base.MPFR.ROUNDING_MODE[])
            return b2
        else
            return BigFloat(b, precision=precision)
        end
    else
        return b
    end
end

function Float128(x::BigFloat)
    if !isfinite(x) || iszero(x)
        return Float128(Float64(x))
    end

    y,k = frexp(x)
    if k >= -16381
        prec = 113
    elseif k >= -16381-112
        prec = 113 + (k + 16381)
    elseif k == -16381-113 && abs(y) > 0.5
        z = reinterepret(Float128, UInt128(1))
    else
        z = reinterepret(Float128, UInt128(0))
    end

    @static if VERSION < v"1.1"
        y = setprecision(BigFloat, prec) do
            BigFloat()
        end
        ccall((:mpfr_set, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32),
                  y, x, Base.MPFR.ROUNDING_MODE[])
    else
        y = BigFloat(x, precision=prec)
    end

    u = zero(UInt128)
    i = cld(prec, sizeof(Culong)*8)
    j = 113
    while i > 0
        j -= sizeof(Culong)*8
        u |= (unsafe_load(y.d, i) % UInt128) << j
        i -= 1
    end
    u &= significand_mask(Float128)
    u |= exponent_half(Float128)
    z = ldexp(reinterpret(Float128, u), y.exp)
    return copysign(z,x)
end

function BigInt(x::Float128)
    !isinteger(x) && throw(InexactError(BigInt, x))
    BigInt(BigFloat(x, precision=precision(Float128)))
end
function Float128(x::BigInt)
    @static if VERSION < v"1.1"
        y = setprecision(BigFloat, precision(Float128)) do
            BigFloat(x)
        end
    else
        y = BigFloat(x, precision=precision(Float128))
    end
    Float128(y)
end

promote_rule(::Type{Float128}, ::Type{Float16}) = Float128
promote_rule(::Type{Float128}, ::Type{Float32}) = Float128
promote_rule(::Type{Float128}, ::Type{Float64}) = Float128
promote_rule(::Type{Float128}, ::Type{<:Integer}) = Float128

#widen(::Type{Float64}) = Float128
widen(::Type{Float128}) = BigFloat

# TODO: need to do this better
function parse(::Type{Float128}, s::AbstractString)
    Float128(ccall((:strtoflt128, libquadmath), Cfloat128, (Cstring, Ptr{Ptr{Cchar}}), s, C_NULL))
end

function string(x::Float128)
    lng = 64
    buf = Array{UInt8}(undef, lng + 1)
    lng = ccall((:quadmath_snprintf,libquadmath), Cint, (Ptr{UInt8}, Csize_t, Ptr{UInt8}, Cfloat128...), buf, lng + 1, "%.35Qe", x)
    return String(resize!(buf, lng))
end

print(io::IO, b::Float128) = print(io, string(b))
show(io::IO, b::Float128) = print(io, string(b))

end # module Quadmath
