using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Petalo

using Test

@testset "Synchronize" begin
    a = (1,2,3,4,  6,7,  9)
    b = (1,2,  4,5,6,  8,9)
    x = (1,2,  4,  6,    9)
    synced = collect(synchronize((a,b), (identity, identity)))
    expected = collect(((n, [n,n]) for n in x))
    @test synced == expected
end

# TODO should be empty when any input is empty
