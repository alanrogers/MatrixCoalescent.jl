using Test
include("../src/LTriMatrix.jl")


@testset "LTriMatrix.jl" begin
    n = 4
    nentries = (n*(n+1)) ÷ 2
    m = LTriMatrix(n, BigInt(1))
    @test size(m) == (n,n)
    @test length(m) == nentries

    for i in 1:n
        for j in 1:i
            m[i,j] = i+j
        end
    end

    @test m[1,1] == 2
    @test m[2,1] == 3
    @test m[3,2] == 5
    @test m[4,4] == 8
    @test_throws BoundsError m[1,2]
end

@testset "BitsLTriMatrix.jl" begin
    n = 4
    nentries = (n*(n+1)) ÷ 2
    m = BitsLTriMatrix(n, 1)
    @test size(m) == (n,n)
    @test length(m) == nentries

    for i in 1:n
        for j in 1:i
            m[i,j] = i+j
        end
    end

    @test m[1,1] == 2
    @test m[2,1] == 3
    @test m[3,2] == 5
    @test m[4,4] == 8
    @test_throws BoundsError m[1,2]
end

