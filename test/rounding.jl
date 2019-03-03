
@testset "rounding" begin
    half128 = Float128(0.5)
    one128 = Float128(1.0)
    a = prevfloat(half128)
    b = half128
    c = (one128 - prevfloat(one128))/2
    d = prevfloat(one128)

    @testset "Default rounding direction, RoundNearest" begin
        @test a + b === one128
        @test - a - b === -one128
        @test a - b === -c
        @test b - a === c
    end
end

@testset "round vs trunc vs floor vs ceil" begin
    x = parse(Float128,"123.456")
    @test round(x, digits=1) ≈ 123.5
    @test round(-x, digits=1) ≈ -123.5
    @test trunc(x, digits=1) ≈ 123.4
    @test trunc(-x, digits=1) ≈ -123.4
    @test ceil(x, digits=1) ≈ 123.5
    @test ceil(-x, digits=1) ≈ -123.4
    @test floor(x, digits=1) ≈ 123.4
    @test floor(-x, digits=1) ≈ -123.5
end

