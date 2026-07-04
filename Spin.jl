module Spin
using .Structs:UnitSpinor
export UnitSpinor_project
UnitSpinor_project(s::UnitSpinor) = s.β / s.α
end

