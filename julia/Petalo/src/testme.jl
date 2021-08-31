using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Petalo

using Test

@testset "trivial stuff" begin
    @test Petalo.double(2) == 4
    @test Petalo.double(π) ≈ 3π rtol = 0.35
end

@testset "Resumable functions" begin
    @test collect(Petalo.onetwothree()) == [1, 2, 3]
end
