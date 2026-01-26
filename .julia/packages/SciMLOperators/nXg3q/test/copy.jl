using SciMLOperators
using LinearAlgebra
using Test

@testset "Copy methods for SciMLOperators" begin
    # Test MatrixOperator
    @testset "MatrixOperator" begin
        A = rand(5, 5)
        L = MatrixOperator(A)
        L_copy = copy(L)
        
        # Modify original
        L.A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.A[1, 1] != 999.0
        @test L_copy.A != L.A
    end
    
    # Test DiagonalOperator (which is a MatrixOperator with Diagonal matrix)
    @testset "DiagonalOperator" begin
        d = rand(5)
        L = DiagonalOperator(d)
        L_copy = copy(L)
        
        # Modify original
        L.A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.A[1, 1] != 999.0
        @test L_copy.A != L.A
    end
    
    # Test ScalarOperator
    @testset "ScalarOperator" begin
        L = ScalarOperator(2.0)
        L_copy = copy(L)
        
        # Modify original
        L.val = 999.0
        
        # Check that copy is not affected
        @test L_copy.val == 2.0
    end
    
    # Test AffineOperator
    @testset "AffineOperator" begin
        A = MatrixOperator(rand(5, 5))
        B = MatrixOperator(rand(5, 5))
        b = rand(5)
        L = AffineOperator(A, B, b)
        L_copy = copy(L)
        
        # Modify original
        L.b[1] = 999.0
        L.A.A[1, 1] = 888.0
        L.B.A[1, 1] = 777.0
        
        # Check that copy is not affected
        @test L_copy.b[1] != 999.0
        @test L_copy.A.A[1, 1] != 888.0
        @test L_copy.B.A[1, 1] != 777.0
    end
    
    # Test ComposedOperator
    @testset "ComposedOperator" begin
        A = MatrixOperator(rand(5, 5))
        B = MatrixOperator(rand(5, 5))
        L = A ∘ B
        L_copy = copy(L)
        
        # Modify original
        L.ops[1].A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.ops[1].A[1, 1] != 999.0
    end
    
    # Test InvertedOperator
    @testset "InvertedOperator" begin
        A = MatrixOperator(rand(5, 5) + 5I)  # Make sure it's invertible
        L = inv(A)
        L_copy = copy(L)
        
        # Modify original
        L.L.A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.L.A[1, 1] != 999.0
    end
    
    # Test TensorProductOperator
    @testset "TensorProductOperator" begin
        A = MatrixOperator(rand(3, 3))
        B = MatrixOperator(rand(2, 2))
        L = kron(A, B)  # Use kron instead of ⊗
        L_copy = copy(L)
        
        # Modify original
        L.ops[1].A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.ops[1].A[1, 1] != 999.0
    end
    
    # Test AdjointOperator  
    @testset "AdjointOperator" begin
        A = MatrixOperator(rand(5, 5))
        L = SciMLOperators.AdjointOperator(A)  # Create AdjointOperator explicitly
        L_copy = copy(L)
        
        # Modify original
        L.L.A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.L.A[1, 1] != 999.0
    end
    
    # Test TransposedOperator
    @testset "TransposedOperator" begin
        A = MatrixOperator(rand(5, 5))
        L = SciMLOperators.TransposedOperator(A)  # Create TransposedOperator explicitly
        L_copy = copy(L)
        
        # Modify original
        L.L.A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.L.A[1, 1] != 999.0
    end
    
    # Test AddedScalarOperator
    @testset "AddedScalarOperator" begin
        α = ScalarOperator(2.0)
        β = ScalarOperator(3.0)
        L = α + β
        L_copy = copy(L)
        
        # Modify original
        L.ops[1].val = 999.0
        
        # Check that copy is not affected
        @test L_copy.ops[1].val == 2.0
    end
    
    # Test ComposedScalarOperator
    @testset "ComposedScalarOperator" begin
        α = ScalarOperator(2.0)
        β = ScalarOperator(3.0)
        L = α * β
        L_copy = copy(L)
        
        # Modify original
        L.ops[1].val = 999.0
        
        # Check that copy is not affected
        @test L_copy.ops[1].val == 2.0
    end
    
    # Test InvertedScalarOperator
    @testset "InvertedScalarOperator" begin
        α = ScalarOperator(2.0)
        L = inv(α)
        L_copy = copy(L)
        
        # Modify original
        L.λ.val = 999.0
        
        # Check that copy is not affected
        @test L_copy.λ.val == 2.0
    end
    
    # Test IdentityOperator (should return self)
    @testset "IdentityOperator" begin
        L = IdentityOperator(5)
        L_copy = copy(L)
        
        # Should be the same object since it's immutable
        @test L === L_copy
    end
    
    # Test NullOperator (should return self)
    @testset "NullOperator" begin
        L = NullOperator(5)
        L_copy = copy(L)
        
        # Should be the same object since it's immutable
        @test L === L_copy
    end
    
    # Test InvertibleOperator
    @testset "InvertibleOperator" begin
        A = rand(5, 5) + 5I  # Make sure it's invertible
        M = MatrixOperator(A)
        F = lu(A)
        L = InvertibleOperator(M, F)
        L_copy = copy(L)
        
        # Modify original
        L.L.A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.L.A[1, 1] != 999.0
    end
    
    # Test ScaledOperator
    @testset "ScaledOperator" begin
        α = ScalarOperator(2.0)
        A = MatrixOperator(rand(5, 5))
        L = α * A
        L_copy = copy(L)
        
        # Modify original
        L.λ.val = 999.0
        L.L.A[1, 1] = 888.0
        
        # Check that copy is not affected
        @test L_copy.λ.val == 2.0
        @test L_copy.L.A[1, 1] != 888.0
    end
    
    # Test AddedOperator
    @testset "AddedOperator" begin
        A = MatrixOperator(rand(5, 5))
        B = MatrixOperator(rand(5, 5))
        L = A + B
        L_copy = copy(L)
        
        # Modify original
        L.ops[1].A[1, 1] = 999.0
        
        # Check that copy is not affected
        @test L_copy.ops[1].A[1, 1] != 999.0
    end
    
    # Test that operators still work correctly after copying
    @testset "Functionality after copy" begin
        # MatrixOperator
        A = rand(5, 5)
        L = MatrixOperator(A)
        L_copy = copy(L)
        v = rand(5)
        @test L * v ≈ L_copy * v
        
        # ScalarOperator
        α = ScalarOperator(2.0)
        α_copy = copy(α)
        @test α * v ≈ α_copy * v
        
        # ComposedOperator
        B = MatrixOperator(rand(5, 5))
        comp = L ∘ B
        comp_copy = copy(comp)
        @test comp * v ≈ comp_copy * v
        
        # AffineOperator
        b = rand(5)
        aff = AffineOperator(L, B, b)
        aff_copy = copy(aff)
        @test aff * v ≈ aff_copy * v
    end
end