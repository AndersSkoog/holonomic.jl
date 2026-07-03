module
using StaticArrays


# complex types
struct C2 <: FieldVector{2,Complex64}
  z1::ComplexF64
  z2::ComplexF64
end

struct CP1Atlas :< FieldVector{2,Complex64}
  ζ::ComplexF64
  ξ::ComplexF64
end

struct UnitSpinor :< FieldVector{2,ComplexF64}
  α::ComplexF64
  β::ComplexF64
end

struct S1 :< FieldVector{2,Float64}
  x::Float64
  z::Float64
end

struct S2 :< FieldVector{3,Float64}
  x::Float64
  y::Float64
  z::Float64
end

struct R3 :< FieldVector{3,Float64}
  x::Float64
  y::Float64
  z::Float64
end

SO3 = SMatrix{3,3,Float64,9}
SU2 = SMatrix{2,2,ComplexF64,4}
su2 = SMatrix{2,2,ComplexF64,4}

struct S3 :< FieldVector{2,Vector{C2}}
  fiber1::Vector{C2}
  fiber2::Vector{C2}
end

#composite types
struct Holomorph2
  index::Int64
  refpoint::ComplexF64
  torangle::Float64
  result::Vector{ComplexF64}
end

Base.:⊗(U::SU2, p::C2) = begin
  v = U * [p.z1; p.z2]
  C2(v[1], v[2])
end

Base.:⊗(U::SU2, fiber::Vector{C2}) = map(p -> U ⊗ p, fiber)








