module Structs
export UnitQuaternion CP1Atlas, R3, C2, S3, UnitSpinor, S2, SO3, SU2, su2, Holomorph2, RollConnection

struct UnitQuaternion <: FieldVector{4,Float64}
  w::Float64
  i::Float64
  j::Float64
  k::Float64
end

struct CP1Atlas :< FieldVector{2,ComplexF64}
  ζ::ComplexF64
  ξ::ComplexF64
end

struct R3 :< FieldVector{3,Float64}
  x::Float64
  y::Float64
  z::Float64
end

struct C2 :< FieldVector{2,Complex64}
  z1::ComplexF64
  z2::ComplexF64
end

struct S3 :< FieldVector{2,Vector{C2}}
  fiber1::Vector{C2}
  fiber2::Vector{C2}
end

struct UnitSpinor :< FieldVector{2,ComplexF64}
  α::ComplexF64
  β::ComplexF64
end

struct S2 :< FieldVector{3,Float64}
  x::Float64
  y::Float64
  z::Float64
end

const SO3 = SMatrix{3,3,Float64,9}
const SU2 = SMatrix{2,2,ComplexF64,4}
const su2 = SMatrix{2,2,ComplexF64,4}

Base.:*(U::SU2, p::C2) = begin
    v = U * [p.z1; p.z2]
    C2(v[1], v[2])
end

struct Holomorph2
  index::Int64
  refpoint::ComplexF64
  torangle::Float64
  frame::SU2
  result::Vector{ComplexF64}
end

struct RollConnection
  contacts::Vector{ComplexF64}
  framesSU2::Vector{SU2}
  framesSO3::Vector{SO3}
  torangles::Vector{Float64}
end

end
