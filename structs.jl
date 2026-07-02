# complex types
struct C2
  z1::ComplexF64
  z2::ComplexF64
end

struct CP1
  z::complexF64
end

struct CP1Atlas
  ζ::ComplexF64
  ξ::ComplexF64
end

struct UnitSpinor
  α::ComplexF64
  β::ComplexF64
end

#cartesian types
struct S1
  x::Float64
  z::Float64
end

struct S2
  x::Float64
  y::Float64
  z::Float64
end

struct R3
  x::Float64
  y::Float64
  z::Float64
end

#matrix types
struct SU2
  U::Matrix{ComplexF64}
end

struct SO3
  R::Matrix{Float64}
end

#array types
struct S3
  fiber1::Array{C2}
  fiber2::Array{C2}
end

#composite types
struct Holomorph2
  index::Int64
  refpoint::ComplexF64
  torangle::Float64
  result::Array{ComplexF64}
end


