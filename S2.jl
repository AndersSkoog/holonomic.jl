module S2
using Random
using .Structs:S2,CP1Atlas
export theta,phi,S2_from_angles,random_closed_sphere_curve,StereoProj,InvStereoProj

function theta(v::S2)
    return acos(clamp(v.z, -1.0, 1.0))
end

function phi(v::S2)
    return atan2(v.y, v.x)
end

function S2_from_angles(θ::Float64,φ::Float64)
    x=cos(φ)*cos(θ)
    y=cos(φ)*sin(θ)
    z=sin(φ)
    n = sqrt(x^2 + y^2 + z^2)
    return S2(x/n,y/n,z/n)
end

function random_closed_sphere_curve(n::Int64=360, k::Int64=5)
    t = range(θ, 2π, length=n+1)[1:end-1]
    θ = zeros(Float64, n)
    φ = zeros(Float64, n)
    for i in 1:k

        Aθ, Bθ = randn(), randn()
        Aφ, Bφ = randn(), randn()

        phaseθ = 2π * rand()
        phaseφ = 2π * rand()

        θ .+= Aθ .* cos.(i .* t .+ phaseθ) .+
            Bθ .* sin.(i .* t .+ phaseθ)

        φ .+= Aφ .* cos.(i .* t .+ phaseφ) .+
            Bφ .* sin.(i .* t .+ phaseφ)
    end

    [S2_from_angles(θ[j],φ[j]) for j in eachindex(θ)]
end

function StereoProj(v::S2)
  ζ = (v.x+v.y*im) / (1-v.z)
  ξ = (v.x-v.y*im) / (1+v.z)
  return CP1Atlas(ζ,ξ)
end

function InvStereoProj(v::ComplexF64)
  m = abs2(v)
  x = (2*v.real)/(1+m)
  y = (2*v.imag)/(1+m)
  z = (1-m)/(1+m)
  n = sqrt(x^2 + y^2 + z^2)
  return S2(x/n,y/n,z/n)
end


end


