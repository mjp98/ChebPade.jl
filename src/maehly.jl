# Based on https://github.com/chebfun/chebfun/blob/master/%40chebfun/chebpade.m

"""
    maehly(F::Fun{<:Chebyshev}, m::Integer, n::Integer)

computes the (m,n) Chebyshev-Pade approximant to F using the Maehly method
"""
function maehly(F::Fun{<:Chebyshev,T}, m::Integer, n::Integer) where T

    # Tolerance for determining negligible singular values.
    tol = 1e-10

    # Get the Chebyshev coefficients and pad if necessary.
    a = extendrandn(copy(coefficients(F)), m + 2n + 1)

   # Form matrix for computing denominator coefficients.
    r = repeat((1:n)', n, 1)
    c = repeat((m+1:m+n), 1, n)

    D = a[c+r.+1] + a[abs.(c - r).+1]

    if n > m
        D += a[1] * diagm(n,n,m=>ones(n - m))
    end
    # D is Hankel + Toeplitz + Diagonal
    # Need to correctly loop a for abs.(c-r)
    # D = Matrix{T}(undef,n,n)
    # @views begin
    #     D  .= Hankel(a[2:n+1],a[m+1:m+n])
    #     D .+= Toeplitz(a[1:n],a[1:n])
    #     if n>m
    #         D .+= a[1] * diagm(n,n,m=>ones(n - m))
    #     end
    # end

    # Solve system for denominator coefficients.
    if (rank(D, tol) < min(size(D)...))
        # If system matrix is singular, reduce degrees first and try again.
        if m > 1
            @info "Singular matrix encountered.\n Reducing m and computing [$(m-1)/$n] approximant."
            return maehly(F, m - 1, n)
        end
        if n > 1
            @info "Singular matrix encountered.\n Reducing n and computing [$m/$(n-1)] approximant"
            return maehly(F, m, n - 1)
        end
        @error "Singular matrix encountered.\n Cannot compute [1/1]  approximation."
        return nothing
    else
        # Otherwise, solve for the denominator coefficients.
        qk = [1; -D \ (2 * a[m+2:m+n+1])]
    end
    # Compute numerator coefficients.
    c = repeat((1:m), 1, n)
    r = repeat((1:n)', m, 1)
    B = a[c+r.+1] + a[abs.(c - r).+1]
    mask = 1:(m+1):min(m, n)*(m+1)
    B[mask] .+= a[1]

    B = [transpose(a[2:n+1]); B]

    # # B is Hankel + Toeplitz + Diagonal
    # Need to correctly loop a for abs.(c-r)
    # B = Matrix{T}(undef,m+1,n)
    # @views begin
    #     B[1,:] .= a[2:n+1]
    #     B[2:m+1,:]  .= Hankel(a[2:m+1],a[m+1:m+n])
    #     B[2:m+1,:] .+= Toeplitz(a[1:m],a[1:n])
    #     B[2:m+1,:] .+= diagm(a[1]*ones(m))
    # end

    pk = qk[1] * a[1:m+1]


    if !isempty(B)
        pk += B * qk[2:n+1]/2
    end

    # Form the outputs.
    p = Fun(space(F), pk)
    q = Fun(space(F), qk)
    return p,q
end
