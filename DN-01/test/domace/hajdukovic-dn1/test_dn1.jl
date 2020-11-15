using Test
using NumMat


@testset "Data types" begin
    #getindex(PM,i,j)
    @testset "getindex" begin
        A = PasovnaMatrika([1,1,1,1],[[3,3,3],[2,2]],[[5,5,5],[6,6]])
        @test A[1,1] == 1
        @test A[2,1] == 5
        @test A[4,1] == 0
        @test A[2,4] == 2
        @test getindex(A,4,4) == 1
        @test getindex(A,2,3) == 3
        @test_throws BoundsError A[5,1]
        @test_throws BoundsError A[0,0]

        B = SpodnjePasovnaMatrika([4,4,4,4],[[-1,3,2],[13.0,2]])
        @test B[1,1] == 4
        @test B[2,1] == -1
        @test B[4,1] == 0
        @test getindex(B,4,4) == 4
        @test getindex(B,3,1) == 13
        @test_throws BoundsError B[5,1]
        @test_throws BoundsError B[0,-1]
        @test_throws BoundsError B[2,4]

        C = ZgornjePasovnaMatrika([4,4,4,4],[[-1,3,2],[13.0,2]])   
        @test C[1,1] == 4
        @test C[2,3] == 3
        @test C[1,4] == 0
        @test getindex(C,4,4) == 4
        @test getindex(C,3,4) == 2
        @test_throws BoundsError C[5,1]
        @test_throws BoundsError C[0,-1]
        @test_throws BoundsError C[4,2]     
end
    @testset "setindex!" begin
        A = PasovnaMatrika([1,1,1,1],[[3,3,3],[2,2]],[[5,5,5],[6,6]])
        A[1,1] = 3
        @test A[1,1] == 3

        A[2,1] = 7
        @test A[2,1] == 7

        setindex!(A,10,3,4)
        @test A[3,4] == 10
        @test_throws BoundsError A[10,10]=5
        @test_throws BoundsError A[0,-1]=10

        B = SpodnjePasovnaMatrika([4,4,4,4],[[-1,3,2],[13.0,2]])
        B[1,1] = -100
        @test B[1,1] == -100
        setindex!(B,10,3,2)
        @test B[3,2] == 10
        @test_throws BoundsError B[1,4] = 1
        @test_throws BoundsError B[-1,-1] = 1

        C = ZgornjePasovnaMatrika([4,4,4,4],[[-1,3,2],[13.0,2]]) 
        C[1,1] = -100
        @test C[1,1] == -100
        setindex!(C,10,2,3)
        @test C[2,3] == 10
        @test_throws BoundsError C[4,1] = 1
        @test_throws BoundsError C[-1,-1] = 1
       
end
end

@testset "Matrix-vector multiplication" begin
    """
    Tests basic functionality of the matrix-vector multiplication. 
    """
    A = PasovnaMatrika(ones(5),[],[])
    B = PasovnaMatrika([1,5,9],[[2,6],[3]],[[4,8],[7]])
    @test A * [1,2,3,4,5] == [1,2,3,4,5]
    @test B * [2,1,3] == [13,31,49]

    """
    Test compares speed(time elapsed) of the matrix-vector multiplication of 
     Regular sparse, n x n matrix and 
     Banded matrix in which only the non-zero elements are stored in the diagonals.
    """

    #Banded matrix representation of regular,sparse n x n matrix with zeros contained  
    # Sparse matrix is represented as PasovnaMatrika in order to make comparisons between same implementation of algorithms. 
    C = NumMat.napolni_vzorcno_matriko!(NumMat.pasovna_nicle(2500,2499,2499),50)
    D = NumMat.napolni_vzorcno_matriko!(NumMat.pasovna_nicle(2500,50,50),50)
    @test ((@elapsed C * ones(2500)) > (@elapsed D * ones(2500))) == true

end


#tests also the  LU decomposition 
@testset "Solving the system" begin
    """
    Tests basic functionality of the matrix-vector multiplication. 
    """
    A = PasovnaMatrika([3,-3,4],[[-2,2]],[[1,2]])
    x = A \ [1.0,1.0,1.0]
    @test isapprox(A*x,[1.0,1.0,1.0],atol=1e-8)
    y = A \ [-5.0,100000.0,-33232.03232]
    @test isapprox(A*y, [-5.0,100000.0,-33232.03232], atol=1e-8)

    B = PasovnaMatrika([3,-3,4],[[-2,10],[1]],[[1,2],[-1]])
    @test_throws ArgumentError("Matrix is not diagonally dominant") B\[1.0,1.0,1.0]

    """
    Test compares speed(time elapsed) of the division from the left of 
     Regular sparse, n x n matrix and 
     Banded matrix in which only the non-zero elements are stored in the diagonals.
    """
    C = NumMat.napolni_vzorcno_matriko!(NumMat.pasovna_nicle(100,99,99),10)
    D = NumMat.napolni_vzorcno_matriko!(NumMat.pasovna_nicle(100,10,10),10)
    @test ((@elapsed C \ ones(100)) > (@elapsed D \ ones(100))) == true

end

