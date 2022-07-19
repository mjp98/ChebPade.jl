const epsreal = eps ∘ float ∘ real

"""
    extendrandn(x::AbstractVector{T},n::Integer)

returns a vector of length n, where the trailing entries are given by randn scaled by eps(float(real(T)))
"""
function extendrandn(x::AbstractVector{T}, n::Integer) where {T}
    length(x) < n ? vcat(x, epsreal(T) * randn(n - length(x))) : x[1:n]
end
