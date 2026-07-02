module Spin
using .CP1
using .S2
export UnitSpinor, project


struct UnitSpinor
    α::ComplexF64
    β::ComplexF64
end

function UnitSpinor(v::S2)
    θ = theta(v)
    φ = phi(v)
    α = cos(θ/2)
    β = exp(im * φ) * sin(θ/2)
    return UnitSpinor(ComplexF64(α), ComplexF64(β))
end


function project(s::UnitSpinor)
    CP1(s.β / s.α)
end

end
