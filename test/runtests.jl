using ChebPade, ApproxFun, LinearAlgebra
using Test

# In part based on https://github.com/chebfun/chebfun/blob/master/tests/chebfun/test_chebpade.m

struct RationalFun{T} <: Function
    p::T
    q::T
end
(r::RationalFun)(z) = r.p(z)/r.q(z)

chebpaderat(args...;kwargs...) = RationalFun(chebpade(args...;kwargs...)...)

@testset "ChebPade.jl" begin
    @testset "Clenshaw-Lord" begin
        tol = 1e-8

        F = Fun(cis,-1..1)
        x = points(space(F),ncoefficients(F))

        P8_8 = chebpaderat(F,8,8)
        P7_4 = chebpaderat(F,7,4)
        P7_9 = chebpaderat(F,7,9)
        P20_1 = chebpaderat(F,20,1)

        Fx = F.(x)

        @test P8_8.(x) ≈ Fx atol = tol
        @test P7_4.(x) ≈ Fx atol = tol
        @test P7_9.(x) ≈ Fx atol = tol
        @test P20_1.(x) ≈ Fx atol = tol


        P8_8 = chebpaderat(F,8,8;method=:maehly)
        P7_4 = chebpaderat(F,7,4;method=:maehly)
        P7_9 = chebpaderat(F,7,9;method=:maehly)
        P20_1 = chebpaderat(F,20,1;method=:maehly)

        Fx = F.(x)

        @test P8_8.(x) ≈ Fx atol = tol
        @test P7_4.(x) ≈ Fx atol = tol
        @test P7_9.(x) ≈ Fx atol = tol
        @test P20_1.(x) ≈ Fx atol = tol
    end

    @testset "adapted from Chebfun" begin
        #% An example which doesn't require degree-reduction.
        dom = -1..3;

        P = Fun(dom,[ 0.5045; -1.3813; 2.1122; 0.0558; -0.6817]);
        Q = Fun(dom,[1; 0.1155; -0.8573; -0.2679; 0.5246 ]);
        R = P./Q
        p,q = chebpade(R, 4, 4;method=:maehly);
        @test norm(P - p) + norm(Q - q) < 100*eps();
        #An example which requires degree-reduction.
        p,q = chebpade(R, 6, 5; method = :maehly);
        @test norm(P - p) + norm(Q - q) < 100*eps();

        # An example by Geddes.
        a = [-464/6375 ; -742/6375 ; 349/12750 ; 512/6375 ; 13/3400 ; 2129/51000 ; 1333/8160 ; 9703/34000];
        b = [-32/85 ; -28/85 ; 1];
        f = Fun(x->evalpoly(x,reverse(a))/evalpoly(x,reverse(b)),-1..1);

        p,q = chebpade(f, 7, 2);
        @test norm(f - p./q, Inf) < 2e-15;

        cp = coefficients(p)
        @test first(cp) ≈ 17/46  atol =  1e-13;

        # Try non-rational and complex-valued functions.
        f = Fun(exp,-1..1);
        p,q = chebpade(f, 2, 3);
        p2,q2 = chebpade(im*f, 2, 3);

        r = RationalFun(p,q)
        r2 = RationalFun(p2,q2)

        @test p2(.3)/p(.3) ≈ im atol=1e-13
        @test r2(-.4)/r(-.4) ≈ im atol=1e-13
    end
end
