#  Contains code that is based in part on Chebfun v5's chebfun/@chebfun/chebpade.m
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

"""
    maehly(F::Fun{<:Chebyshev}, m::Integer, n::Integer)

computes the (m,n) Chebyshev-Pade approximant to F using the Maehly method
"""
function maehly(
    F::Fun{<:Chebyshev,T},
    m::Integer,
    n::Integer;
    # Tolerance for determining negligible singular values.
    tol= 10000*epsreal(T)
    ) where T

    pk,qk = maehly(coefficients(F),m,n;tol)

    # Form the outputs.
    p = Fun(space(F), pk)
    q = Fun(space(F), qk)
    return p,q
end

function maehly(
    F::AbstractVector{T},
    m::Integer,
    n::Integer;
    # Tolerance for determining negligible singular values.
    tol= 10000*epsreal(T)
    ) where T
    # Pad coefficients with noise if necessary.
    a = padwithnoise(F, m + 2n + 1)
    # Form matrix for computing denominator coefficients.
    D = [maehly_matrix(a,i,j) for i in m+1:m+n, j in 1:n]
    # Solve system for denominator coefficients.
    if (rank(D, tol) < min(size(D)...))
        # If system matrix is singular, reduce degrees first and try again.
        m > 1 && return maehly(F, m - 1, n)
        n > 1 && return maehly(F, m, n - 1)
        @error "Singular matrix encountered.\n Cannot compute [1/1]  approximation."
        return nothing
    else
        # Otherwise, solve for the denominator coefficients.
       qk = [1; -D \ (2 * a[m+2:m+n+1])]
    end
    # Compute numerator coefficients.
    B = [maehly_matrix(a,i,j) for i in 1:m, j in 1:n]
    B = [transpose(a[2:n+1]); B]
    pk = qk[1] * a[1:m+1]
    if !isempty(B)
        pk += B * qk[2:n+1]/2
    end
    return pk, qk
end

maehly_matrix(x,i,j) = x[i+j+1] + x[abs(i-j)+1] + (i==j)*x[1]
