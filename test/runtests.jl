using MatrixCoalescent
using Test

@testset "Rational" begin
    n=3
    mc = MatCoal(Rational{Int}, n)
    @test mc.nLineages == 3
    @test dim(mc) == n-1
    @test mc.beta[1] == 1//1
    @test mc.beta[2] == 3//1

    @test mc.rbeta[1] == 1//1
    @test mc.rbeta[2] == 1//3

    @test mc.gmat[1,1] == 3//2
    @test mc.gmat[1,2] == -3//2
    @test mc.gmat[2,1] == 0//1
    @test mc.gmat[2,2] == 1//1

    @test mc.hmat[1,1] == -3//2
    @test mc.hmat[1,2] == 1//2
    @test mc.hmat[2,1] == 0//1
    @test mc.hmat[2,2] == -1//3
end

@testset "Float64" begin
    n=3
    mc = MatCoal(Float64, n)
    @test mc.nLineages == 3
    @test dim(mc) == n-1
    @test mc.beta[1] == 1.0
    @test mc.beta[2] == 3.0

    @test mc.rbeta[1] == 1.0
    @test mc.rbeta[2] == 1.0/3.0

    @test mc.gmat[1,1] == 1.5
    @test mc.gmat[1,2] == -3.0/2.0
    @test mc.gmat[2,1] == 0.0/1.0
    @test mc.gmat[2,2] == 1.0

    @test mc.hmat[1,1] == -1.5
    @test mc.hmat[1,2] == 0.5
    @test mc.hmat[2,1] == 0.0
    @test mc.hmat[2,2] == -1.0/3.0

    eig = Vector{Float64}(undef, n-1)
    v = 1.0
    eigenvals!(eig, v, mc)

    @test eig[1] == exp(-1)
    @test eig[2] == exp(-3)

    pr = Vector{Float64}(undef, n-1)
    project!(pr, eig, mc)

    @test pr[1] ≈ 1.5 * (exp(-1) - exp(-3))
    @test pr[2] ≈ exp(-3)

    elen = Vector{Float64}(undef, n-1)
    interval_lengths!(elen, eig, mc)

    elen[1] ≈ 1.0 - 1.5*exp(-1) + 0.5*exp(-3)
    elen[2] ≈ 1.0/3.0 - exp(-3)/3.0
end

"""
For a given number of lineages (nLin) and across a range of time values,
calculate the projected probability vector (using project!) and the
expected interval lengths (using interval_lengths!) using 64-bit
floats and also using 256-bit floats. Return the largest absolute
difference between the low-precision and the high-precision versions
of this calculation.
"""
function estimate_error(nLin::Int)
    mcl = MatCoal(Float64, nLin)  # low precision
    mch = MatCoal(BigFloat, nLin) # high precision

    el = Vector{Float64}(undef, nLin-1)
    eh = Vector{BigFloat}(undef, nLin-1)
    xl = Vector{Float64}(undef, nLin-1)
    xh = Vector{BigFloat}(undef, nLin-1)

    err = 0.0

    for v in 0.0 : 0.5 : 9.5
        eigenvals!(el, v, mcl)
        eigenvals!(eh, BigFloat(v), mch)
        project!(xl, el, mcl)
        project!(xh, eh, mch)
        for (l,h) in zip(xl, xh)
            err = max(err, abs(l - Float64(h)))
        end
        interval_lengths!(xl, el, mcl)
        interval_lengths!(xh, eh, mch)
        for (l,h) in zip(xl, xh)
            err = max(err, abs(l - Float64(h)))
        end
    end
    return err
end

@testset "Accuracy" begin
    @test estimate_error(2)  < eps(Float64)
    @test estimate_error(4)  < eps(Float64)
    @test estimate_error(8)  < 4e-15
    @test estimate_error(16) < 3.5e-13
    @test estimate_error(32) < 4e-8
    @test estimate_error(34) < 6e-8
    @test estimate_error(35) < 3e-7
    @test_throws OverflowError estimate_error(36)
end

