# Based on https://github.com/chebfun/chebfun/blob/master/%40chebfun/chebpade.m

"""
    clenshawlord(F::Fun{<:Chebyshev}, m::Integer, n::Integer)

computes the (m,n) Chebyshev-Pade approximant to F using the Clenshaw-Lord method
"""
function clenshawlord(F::Fun{<:Chebyshev}, m::Integer, n::Integer)

    c = extendrandn(copy(coefficients(F)), m + 2n + 1)

    l = max(m, n)                   # Temporary degree variable in case m < n.

    c[1] *= 2

    # Set up and solve Hankel system for denominator Laurent-Pade coefficients.
    top = c[abs.((m-n+1:m)).+1]     # Top row of Hankel system.
    bot = c[(m:m+n-1).+1]           # Bottom row of Hankel system.
    H = Hankel(top, bot)            # Hankel matrix
    rhs = c[(m+1:m+n).+1]           # RHS of Hankel system.

    # Use convolution to compute numerator Laurent-Pade coefficients.
    β = n == 0 ? [1] : reverse([-H \ rhs; 1])
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

    p = Fun(space(F), pk)
    q = Fun(space(F), qk)
    return p,q
end
