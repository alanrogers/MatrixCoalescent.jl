using Test
include("../src/UTriMatrix.jl")

@testset "UTriMatrix.jl" begin
    n = 4
    nentries = (n*(n+1)) ÷ 2
    m = UTriMatrix(n, BigInt(1))
    @test size(m) == (n,n)
    @test length(m) == nentries

    for i in 1:n
        for j in i:n
            m[i,j] = i+j
        end
    end

    @test m[1,1] == 2
    @test m[1,2] == 3
    @test m[2,3] == 5
    @test m[4,4] == 8
    @test_throws BoundsError m[2,1]
end

@testset "BitsUTriMatrix.jl" begin
    n = 4
    nentries = (n*(n+1)) ÷ 2
    m = BitsUTriMatrix(n, 1)
    @test size(m) == (n,n)
    @test length(m) == nentries

    for i in 1:n
        for j in i:n
            m[i,j] = i+j
        end
    end

    @test m[1,1] == 2
    @test m[1,2] == 3
    @test m[2,3] == 5
    @test m[4,4] == 8
    @test_throws BoundsError m[2,1]
end

