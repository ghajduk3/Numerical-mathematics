using Test
using NumMat

@testset "dn3" begin
    @testset "Chebyshev nodes" begin
        a1,b1,x1,n1 = (0.001,5,2.3,5)
        a2,b2,x2,n2 = (2,15,4.3,10)
        a3,b3,x3,n3 = (10,15,12.4,15)
        cheb1 = cos((x1-1))*pi/n1
        cheb2 = cos((x2-1))*pi/n2
        cheb3 = cos((x3-1))*pi/n3
        @test isapprox(NumMat.cheb_nodes(a1,b1,cheb1),2.9206021401,atol=10e-6)
        @test isapprox(NumMat.cheb_nodes(a2,b2,cheb2),6.48353152560, atol=10e-6)
        @test isapprox(NumMat.cheb_nodes(a3,b3,cheb3),12.70603133, atol=10e-6)

    end

    @testset "Barycentric weights" begin
        @test NumMat.weight(5)[1] == 0.5
        @test NumMat.weight(5)[end] == -0.5
        @test NumMat.weight(5)[2] == -1.0
        @test NumMat.weight(5)[3] == 1.0
    
    end

    @testset "Barycentric interpolation" begin
    # Function to be interpolated x->exp((-x)^2) on interval [-1,1]
        for i=1:1000
            i_pol = NumMat.inter_polynom(x->exp((-x)^2),-1,1,i)
            values = Array{Float64}(undef,50)
            for j=1:50
                values[j] = NumMat.inter_value(i_pol,i_pol.X[j])
            end
            if  sum(abs.(i_pol.Y.-values))/length(values) < 10e-6
                # "Polynomial degree"
                @test i==10
                # Precision 
                @test isapprox(values,map(x->exp((-x)^2),i_pol.X),atol=10e-6)
                break
            end    
        end
        
        # Function to be interpolated x->sin(x)/x on interval [0.001,10]
        for i=1:1000
            i_pol = NumMat.inter_polynom(x->sin(x)/x,0.001,10,i)
            values = Array{Float64}(undef,50)
            for j=1:50
                values[j] = NumMat.inter_value(i_pol,i_pol.X[j])
            end
            if  sum(abs.(i_pol.Y.-values))/length(values) < 10e-6
                # "Polynomial degree"
                @test i==14
                # Precision 
                @test isapprox(values,map(x->sin(x)/x,i_pol.X),atol=10e-6)
                break
            end    
        end 


    end



end