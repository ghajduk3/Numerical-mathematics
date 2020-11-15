using Test
using NumMat


@testset "dn2" begin
    RM =  RazprsenaMatrika([[5,0,2,1,0],[0,5,3,2,0],[0,0,0,0,0],[1,0,6,0,0],[0,5,7,0,0]])
    regular = [5.0 0.0 2.0 1.0 0.0; 0.0 5.0 3.0 2.0 0.0; 0.0 0.0 0.0 0.0 0.0 ; 1.0 0.0 6.0 0.0 0.0 ; 0.0 5.0 7.0 0.0 0.0]
    b = [2.3,1.5,3.0,3.1,0.53]

    regular1,RM1 = NumMat.fill_sample_matrix!([zeros(1000) for x=1:1000], 1000)
    b1 = rand(Int,1000)

    regular2,RM2 = NumMat.fill_sample_matrix!([zeros(10000) for x=1:10000], 10000)
    b2 = rand(Int,10000)

    SOR1 = RazprsenaMatrika([[2.0,-1.0,0.0],[-1.0,2.0,-1.0],[0.0,-1.0,2.0]])
    regular3 = [2.0 -1.0 0.0 ; -1.0 2.0 -1.0 ; 0.0 -1.0 2.0]
    b3 = [0.0,1.0,2.0]

    x = [1.0,2.0,3.0]
    b4 = regular3 * x

    x = [0.0, 0.0 , 0.0]
    b5 = regular3 * x


    @testset "Indexing" begin
        @testset "getindex" begin
            @test RM[1,1] == 5
            @test RM[3,1] == 0
            @test RM[5,3] == 7
            @test getindex(RM,5,3) == RM[5,3]
            @test getindex(RM,5,5) == 0
            @test_throws BoundsError RM[-5,3]
            @test_throws BoundsError RM[10,10]
    end
        @testset "setindex!" begin
            RM[1,1] = 3
            @test RM[1,1] == 3

            RM[2,1] = 7
            @test RM[2,1] == 7

            setindex!(RM,10,3,4)
            @test RM[3,4] == 10
            @test_throws BoundsError RM[10,10]=5
            @test_throws BoundsError RM[0,-1]=10

    end
        @testset "First index" begin
            @test firstindex(RM) == 1
            @test firstindex(RM1) == 1
            @test firstindex(RM2) == 1
    end
        @testset "Last index" begin
            @test lastindex(RM) == RM.n
            @test lastindex(RM1) == RM1.n
            @test lastindex(RM2) == RM2.n
    end
end
    @testset "Matrix-vector multiplication" begin
        @test RM * ones(5) == [6.0,17.0,10.0,7.0,12.0]
        @test RM1 * ones(1000) == regular1 * ones(1000)
        @test @elapsed(RM2 * ones(10000)) < @elapsed(regular2 * ones(10000)) 
        B = [1.5 2 -4; 2 -1 -3; -4 -3 5]
    end   

    @testset "SOR" begin
       @test isapprox(sor(SOR1,b3,[0.0,0.5,1.0],1.2,1e-10)[1] ,regular3 \ b3,atol=10e-8)
       @test isapprox(sor(SOR1,b4,[0.0,0.0,0.0],1.5,1e-10)[1] ,regular3 \ b4,atol=10e-8)
       @test isapprox(sor(SOR1,b5,[0.0,0.0,0.0],1.5,1e-10)[1] ,regular3 \ b5,atol=10e-8) 
    end
end