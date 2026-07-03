module R3
using StaticArrays
export R3

struct R3 :< FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end



end
