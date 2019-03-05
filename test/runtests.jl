using Test
using Quadmath



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
