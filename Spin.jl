module Spin
using .Structs:UnitSpinor
export UnitSpinor_project, UnitSpinor_from_S2

function UnitSpinor_from_S2(v::S2)
    θ = theta(v)
    φ = phi(v)
    α = cos(θ/2)
    β = exp(im * φ) * sin(θ/2)
    return UnitSpinor(ComplexF64(α), ComplexF64(β))
end

UnitSpinor_project(s::UnitSpinor) = s.β / s.α
end

