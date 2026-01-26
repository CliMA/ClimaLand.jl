using PolygonOps
using Test
using StaticArrays

edge(p1, p2) = (SVector(p1),SVector(p2))

@testset "inpolygon" begin
    poly1 = SVector{2,Int}[(0,0),(0,10),(10,10),(10,0),(0,0)]

    @testset "Hao Sun" begin
        algo = HaoSun()
        for poly in (poly1, reverse(poly1))
            @test inpolygon(SVector(5,5),poly, algo) == 1 #in
            @test inpolygon(SVector(-5,-5),poly, algo) == 0 #out
            @test inpolygon(SVector(0,5),poly, algo) == -1 #on
            @test inpolygon(SVector(5,10),poly, algo) == -1 #on
            @test inpolygon(SVector(10,10),poly, algo) == -1 #on
            @test inpolygon(SVector(0,0),poly, algo) == -1 #on

            @test inpolygon(SVector(5,5),poly, algo, in=true,on=true,out=false) #in
            @test !inpolygon(SVector(-5,-5),poly, algo, in=true,on=true,out=false) #out
            @test inpolygon(SVector(0,5),poly, algo, in=true,on=true,out=false) #on
            @test inpolygon(SVector(5,10),poly, algo, in=true,on=true,out=false) #on
            @test inpolygon(SVector(10,10),poly, algo, in=true,on=true,out=false) #on
            @test inpolygon(SVector(0,0),poly, algo, in=true,on=true,out=false) #on
        end
    end
    @testset "Hormann Agathos" begin
        algo = HormannAgathos()

        poly = reverse(poly1)
        @test inpolygon(SVector(5,5),poly, algo) == 1 #in
        @test inpolygon(SVector(-5,-5),poly, algo) == 0 #out
        #@test_broken inpolygon(SVector(0,5),poly, algo) == -1 #on
        @test inpolygon(SVector(5,10),poly, algo) == -1 #on
        @test inpolygon(SVector(10,10),poly, algo) == -1 #on
        @test inpolygon(SVector(0,0),poly, algo) == -1 #on
    end
    @testset "API" begin
        poly = poly1
        @test inpolygon(SVector(5,5),poly, in=true,on=true,out=false) #in
        @test !inpolygon(SVector(-5,-5),poly, in=true,on=true,out=false) #out
        @test inpolygon(SVector(0,5),poly, in=true,on=true,out=false) #on
        @test inpolygon(SVector(5,10),poly, in=true,on=true,out=false) #on
        @test inpolygon(SVector(10,10),poly, in=true,on=true,out=false) #on
        @test inpolygon(SVector(0,0),poly, in=0x0,on=0x0,out=0x1) === 0x0 #on
    end
end

@testset "area/centroid" begin
    poly1 = SVector{2,Int}[(0,0),(0,10),(10,10),(10,0),(0,0)]
    poly2 = poly1.*5

    @testset "area" begin
        @test PolygonOps.area(poly1) == -100
        @test PolygonOps.area(reverse(poly1)) == 100
        @test PolygonOps.area(poly2) == -2500
        @test PolygonOps.area(reverse(poly2)) == 2500
    end

    @testset "centroid" begin
        @test PolygonOps.centroid(poly1) == SVector(5,5)
        @test PolygonOps.centroid(reverse(poly1)) == SVector(5,5)
        @test PolygonOps.centroid(poly2) == SVector(25,25)
        @test PolygonOps.centroid(reverse(poly2)) == SVector(25,25)
    end
end
