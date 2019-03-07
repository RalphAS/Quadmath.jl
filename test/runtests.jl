using Test
using Quadmath

using InteractiveUtils
@show reinterpret(UInt128,Float128(3.0))

@show reinterpret(UInt128,Float128(3.0) + Float128(4.0))
#=
using Quadmath: quadoplib, Cfloat128
function foo(x::Float128,y::Float128)
    Float128(ccall((:__addtf3,quadoplib),
                   Cfloat128, (Ref{Cfloat128},Ref{Cfloat128}), x, y))
end
@show reinterpret(UInt128,foo(Float128(3.0),Float128(3.0)))
=#
let y = reinterpret(Float128, 0x40008000000000000000000000000000)
    @code_native frexp(y)
    @code_native ldexp(y, 2)
    @code_native ==(y,y)
#    @code_native foo(y,y)
end


@testset "fp decomp" begin
    y = Float128(2.0)
    z = ldexp(Float128(0.5), 2)
    @test z == y
    x,n = frexp(y)
    @test x == Float128(0.5)
    @test n == 2
end

@testset "conversion $T" for T in (Float64, Int32, Int64, BigFloat, BigInt)
    @test Float128(T(1)) + Float128(T(2)) == Float128(T(3))
    @test Float128(T(1)) + Float128(T(2)) <= Float128(T(3))
    @test Float128(T(1)) + Float128(T(2)) != Float128(T(4))
    @test Float128(T(1)) + Float128(T(2)) < Float128(T(4))
    if isbitstype(T)
        @test T(Float128(T(1)) + Float128(T(2))) === T(3)
    else
        @test T(Float128(T(1)) + Float128(T(2))) == T(3)
    end
end

@test Base.exponent_one(Float128) == reinterpret(UInt128, Float128(1.0))

@testset "BigFloat" begin
    x = parse(Float128, "0.1")
    y = parse(Float128, "0.2")
    @test Float64(x+y) == Float64(BigFloat(x) + BigFloat(y))
    @test x+y == Float128(BigFloat(x) + BigFloat(y))
end

@testset "BigInt" begin
    x = parse(Float128, "100.0")
    y = parse(Float128, "25.0")
    @test Float64(x+y) == Float64(BigInt(x) + BigInt(y))
    @test x+y == Float128(BigInt(x) + BigInt(y))    
end

@testset "flipsign" begin
    x = Float128( 2.0)
    y = Float128(-2.0)
    @test x == flipsign(y, -one(Float128))
    @test y == flipsign(y,  1)
end

@testset "modf" begin
    x = Float128(pi)
    fpart, ipart = modf(x)
    @test x == ipart + fpart
    @test signbit(fpart) == signbit(ipart) == false

    y = Float128(-pi)
    fpart, ipart = modf(y)
    @test y == ipart + fpart
    @test signbit(fpart) == signbit(ipart) == true
    
    z = x^3
    fpart, ipart = modf(x) .+ modf(y)
    @test x+y == ipart+fpart
end
