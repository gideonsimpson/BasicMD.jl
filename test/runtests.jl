using Test
using BasicMD
using TestLandscapes
using ForwardDiff
using Random
using Statistics

@testset "RWM" begin
    @test include("rwm/rwm1.jl")
    @test include("rwm/rwm2.jl")
    @test include("rwm/rwm3.jl")
    @test include("rwm/rwm4.jl")
    @test include("rwm/rwm5.jl")
end

