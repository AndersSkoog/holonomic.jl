const TNum        = Float64
const TAng        = Float64
const TCplx       = ComplexF64

const TNum2       = Tuple{TNum,TNum}
const TAng2       = Tuple{TAng,TAng}
const TCplx2      = Tuple{TCplx,TCplx}

const TVec2       = SVector{2,TNum}
const TVec3       = SVector{3,TNum}

const TSO3       = SMatrix{3, 3, Float64, 9}
const TSU2       = SMatrix{2, 2, ComplexF64, 4}
const TSE3       = SMatrix{4, 4, Float64, 16}

const TArrNum     = Vector{TNum}
const TArrAng     = Vector{TAng}
const TArrCplx    = Vector{TCplx}
const TArrNum2    = Vector{TNum2}
const TArrAng2    = Vector{TAng2}
const TArrCplx2   = Vector{TCplx2}
const TArrVec2    = Vector{TVec2}
const TArrVec3    = Vector{TVec3}
const TArrSO3     = Vector{TSO3}
const TArrSU2     = Vector{TSU2}
const TArrSE3     = Vector{TSE3}

const TMöbiusCoef        = Tuple{TCplx,TCplx,TCplx,TCplx}
const TFourierSeriesCoef = Tuple{TArrNum,TArrNum}
const THopfFibre         = TArrCplx2
const THopfLink          = Tuple{THopfFibre,THopfFibre}
const TPlane3            = Tuple{TVec3,TVec3}

const TSCurveCoef   = Tuple{TArrNum,TArrNum,TArrNum,TArrNum}

struct TSCurve
    coef::TSCurveCoef
    azi::TArrAng
    pol::TArrAng
    pts::TArrVec3
end

struct TSConnection
    contacts::TSCurve
    form::Vector{Tuple{TSO3,TVec3}}
    pos::TArrVec3
    orient::TArrSO3
    dev::TArrCplx
end




