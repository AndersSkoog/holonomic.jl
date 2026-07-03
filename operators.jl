
Base.:⊗(U::SU2, p::C2) = begin
  v = U * [p.z1; p.z2]
  C2(v[1], v[2])
end

Base.:⊗(U::SU2, fiber::Vector{C2}) = map(p -> U ⊗ p, fiber)

