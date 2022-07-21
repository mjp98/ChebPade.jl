
#  Contains code that is based in part on Chebfun v5's chebfun/chebfun/test_chebpade.m,
# which is distributed with the following license:

# Copyright (c) 2015, The Chancellor, Masters and Scholars of the University
# of Oxford, and the Chebfun Developers. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Oxford nor the names of its
#       contributors may be used to endorse or promote products derived from
#       this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


using ChebPade, ApproxFunOrthogonalPolynomials, LinearAlgebra
using Test

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

    @testset "algorithm specification error" begin
        @test_throws AssertionError chebpade(Fun(zero,-1..1), 4, 4;method=:rand);
    end


    @testset "maehly singular" begin
        # example where degree reduction in denominator is required
        dom = -1..1
        P = Fun(dom,[1]);
        Q = Fun(dom,[2]);
        R = P./Q
        p,q = chebpade(R,1,5;method=:maehly)
        @test norm(R - p./q) < 100eps()
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
