using MyStatsPackage
using Test

@testset "MyStatsPackage.jl" begin
    # Write your tests here.
    @test mean_value([1, 2, 3, 4, 5]) == 3
    @test median_value([5, 2, -4, 4, 1]) == 2
    @test median_value([2, 5, 3, 4]) == 3.5
    @test (variance_value([7,1,6,12,6,2,7,4]) - 11.6964) <= 0.01 
    @test (standard_deviation([7,1,6,12,6,2,7,4]) - 3.42) <= 0.01 

end
