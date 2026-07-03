module CP1
using .S2
export CP1Atlas

struct CP1Atlas :< FieldVector{2,ComplexF64}
    ζ::ComplexF64
    ξ::ComplexF64
end

end
