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
    clenshawlord(F::Fun{<:Chebyshev}, m::Integer, n::Integer)

computes the (m,n) Chebyshev-Pade approximant to F using the Clenshaw-Lord method
"""
function clenshawlord(F::Fun{<:Chebyshev,T}, m::Integer, n::Integer) where T
    pk,qk = clenshawlord(coefficients(F),m,n)
    p = Fun(space(F), pk)
    q = Fun(space(F), qk)
    return p,q
end

function clenshawlord(F::Vector{T}, m::Integer, n::Integer) where T

    c = padwithnoise(F, m + 2n + 1)

    l = max(m, n)                   # Temporary degree variable in case m < n.

    c[1] *= 2

    # Set up and solve Hankel system for denominator Laurent-Pade coefficients.
    top = c[abs.((m-n+1:m)).+1]     # Top row of Hankel system.
    bot = c[(m:m+n-1).+1]           # Bottom row of Hankel system.
    H = Hankel(top, bot)            # Hankel matrix
    rhs = c[(m+1:m+n).+1]           # RHS of Hankel system.

    # Use convolution to compute numerator Laurent-Pade coefficients.
    β = n == 0 ? [one(T)] : reverse([-H \ rhs; 1])
    c[1] /= 2
    α = conv(c[1:l+1], β)
    D = PaddedView(0, [a * b for a in α, b in β], (l + 1, l + 1))

    # numerator
    pk = [sum(diag(D, k)) + sum(diag(D, -k)) for k = 0:m]
    pk[1] /= 2

    # denominator
    qk = @views [dot(β[k:end]', β[1:n+2-k]) for k = 1:n+1]

    # Normalise
    pk ./= qk[1]
    qk .= (2 / qk[1]) * qk
    qk[1] /= 2
    return pk, qk
end
