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

@testset "EM" begin
    @test include("em/em1.jl")
    @test include("em/em2.jl")
    @test include("em/em3.jl")
    @test include("em/em4.jl")
    @test include("em/em5.jl")
end

@testset "MALA" begin
    @test include("mala/mala1.jl")
    @test include("mala/mala2.jl")
    @test include("mala/mala3.jl")
    @test include("mala/mala4.jl")
    @test include("mala/mala5.jl")
end

@testset "HMC" begin
    @test include("hmc/hmc1.jl")
    @test include("hmc/hmc2.jl")
    @test include("hmc/hmc3.jl")
    @test include("hmc/hmc4.jl")
    @test include("hmc/hmc5.jl")
end

@testset "BBK" begin
    @test include("bbk/bbk1.jl")
    @test include("bbk/bbk2.jl")
    @test include("bbk/bbk3.jl")
    @test include("bbk/bbk4.jl")
    @test include("bbk/bbk5.jl")
end

@testset "GJF" begin
    @test include("gjf/gjf1.jl")
    @test include("gjf/gjf2.jl")
    @test include("gjf/gjf3.jl")
    @test include("gjf/gjf4.jl")
    @test include("gjf/gjf5.jl")
end
