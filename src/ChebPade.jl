module ChebPade

using LinearAlgebra, ToeplitzMatrices, PaddedViews
using ApproxFunOrthogonalPolynomials
using DSP

export chebpade

"""
    chebpade(f::Fun{<:Chebyshev},m::Integer,n::Integer;method=:clenshawlord)

computes the (m,n) Chebyshev-Pade approximation to f using method :clenshawlord or :maehly
"""
function chebpade(f::Fun{<:Chebyshev}, m::Integer, n::Integer; method=:clenshawlord)
    @assert method in [:clenshawlord,:maehly] "$method not implemented.\n Set method = :clenshawlord or :maehly"
    method == :clenshawlord && return clenshawlord(f, m, n)
    method == :maehly && return maehly(f, m, n)
end

include("util.jl")
include("clenshawlord.jl")
include("maehly.jl")

end
