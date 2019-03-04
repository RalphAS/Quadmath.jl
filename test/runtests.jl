using Test
using Quadmath

using Quadmath: quadoplib, Cfloat128
@noinline function foobar(y::Float64)
    x = 3.0*y
    r = Base.llvmcall("""
            %f = inttoptr i64 %1 to <2 x double> (double)*
            %vv = call x86_vectorcallcc <2 x double>  %f(double %0)
            ret <2 x double> %vv""", Cfloat128, Tuple{Cdouble,Ptr{Cvoid}},
                 x, cglobal((:__extenddftf2,quadoplib)))
    Float128(r)
end
using InteractiveUtils
@code_native foobar(1.0)
print("foobar(1.0) returns ")
display(reinterpret(UInt128,foobar(1.0)))
println()

function baz(x::Float128, y::Float128)
    r = Base.llvmcall("""%f = inttoptr i64 %2 to <2 x double> (<2 x double>, <2 x double>)*
                    %vv = call x86_vectorcallcc <2 x double> %f(<2 x double> %0, <2 x double> %1)
                    ret <2 x double> %vv""",
                 Cfloat128, Tuple{Cfloat128,Cfloat128,Ptr{Cvoid}},
                 x.data, y.data, cglobal((:__addtf3,quadoplib)))
    Float128(r)
end
x3 = reinterpret(Float128, 0x40008000000000000000000000000000)
x2 = reinterpret(Float128, 0x40000000000000000000000000000000)
@code_native baz(x2,x3)
print("baz(x2,x3) returns ")
display(reinterpret(UInt128,baz(x2,x3)))
println()


@info "starting test run"
    for j in (1.0,2.0,3.0)
        print("Float128($j) is ")
        display(reinterpret(UInt128,Float128(j)))
        println()
    end
    print("sum is ")
    display(reinterpret(UInt128,Float128(1.0) + Float128(2.0)))
    println()
    print("return trip for 2 is ")
    display(reinterpret(UInt64,Float64(Float128(2.0))))
    println()


@testset "conversion $T" for T in (Float64, Float32, Int32, Int64, BigFloat, BigInt)
    @test Float128(T(1)) + Float128(T(2)) == Float128(T(3))
    @test Float128(T(1)) + Float128(T(2)) <= Float128(T(3))
    @test Float128(T(1)) + Float128(T(2)) != Float128(T(4))
    @test Float128(T(1)) + Float128(T(2)) < Float128(T(4))
    for j=1:3
        print("Float128($T($j)) is ")
        display(reinterpret(UInt128,Float128(T(j))))
        println()
    end
    print("sum is ")
    display(reinterpret(UInt128,Float128(T(1)) + Float128(T(2))))
    println()
    if T == Float64
        print("return trip for 2 is ")
        display(reinterpret(UInt64,T(Float128(T(2)))))
        println()
    end
    if isbitstype(T)
        @test T(Float128(T(1)) + Float128(T(2))) === T(3)
    else
        @test T(Float128(T(1)) + Float128(T(2))) == T(3)
    end
end

@testset "conversion Int128" begin
    @test Float128(Int128(123)) == Float128(123.0)
    x = Int128(UInt128(1) << 115)
    @test Float128(x) == Float128(2.0^115)
    @test Float128(-x) == Float128(-2.0^115)
end
@testset "conversion UInt128" begin
    @test Float128(UInt128(123)) == Float128(123.0)
    x = UInt128(1) << 115
    @test Float128(x) == Float128(2.0^115)
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

include("rounding.jl")
