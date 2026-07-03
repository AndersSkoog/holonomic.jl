module Spin
using .CP1:CP1
using .S2:S2,theta,phi
export UnitSpinor, project

struct UnitSpinor :< FieldVector{2,ComplexF64}
    α::ComplexF64
    β::ComplexF64
    function UnitSpinor(v::S2)
        θ = theta(v)
        φ = phi(v)
        α = cos(θ/2)
        β = exp(im * φ) * sin(θ/2)
        new(ComplexF64(α), ComplexF64(β))
    end
end

project(s::UnitSpinor) = s.β / s.α

end
