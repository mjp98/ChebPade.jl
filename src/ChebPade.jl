module ChebPade

using LinearAlgebra, ToeplitzMatrices, PaddedViews
using ApproxFun
using DSP

export chebpade

"""
    chebpade(f::Fun{<:Chebyshev},m::Integer,n::Integer;method=:clenshawlord)

computes the (m,n) Chebyshev-Pade approximation to f using method :clenshawlord or :maehly
"""
function chebpade(f::Fun{<:Chebyshev}, m::Integer, n::Integer; method=:clenshawlord)
    method == :clenshawlord && return clenshawlord(f, m, n)
    method == :maehly && return maehly(f, m, n)
    @error "$method not implemented. Set method = :clenshawlord or :maehly"
end

include("util.jl")
include("clenshawlord.jl")
include("maehly.jl")

end
