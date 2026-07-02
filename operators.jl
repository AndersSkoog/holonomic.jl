
Base.:⊗(U::SU2, p::C2) = begin
  v = U.U * [p.z1; p.z2]
  C2(v[1], v[2])
end

Base.:⊗(U::SU2, fiber::Vector{C2}) = map(p -> U ⊗ p, fiber)


#orient 3d unit vector
Base.:⊗(A::S2,B::SO3) = begin
  arr = A.R * [B.x,B.y,B.z]
  S2(arr[1],arr[2],arr[3])
end

#SO3 matrix multiplication
Base.:⊗(A::SO3,B::SO3) = SO3(A.R*B.R)
#SO3 matrix multiplication
Base.:⊗(A::SU2, B::SU2) = SU2(A.U * B.U)

Base.adjoint(A::SU2) = SU2(adjoint(A.U))
Base.inv(A::SU2) = SU2(adjoint(A.U))  # valid for SU(2)


