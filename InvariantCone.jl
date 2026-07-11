module InvariantCone
using .Structs:S2
export cap_circle
function cap_circle(v::S2)
    n = [0.0,0.0,1.0]
    s = [v.x,v.y,v.z]
    c = (4/5)*s
    r = sin(pi/5)
    θ = theta(v)
    φ = phi(v)
    c1=r*sin(θ)
    c2=r*sin(φ)*cos(θ)
    c3=r*cos(θ)
    c4=r*sin(φ)*sin(θ)
    c5=r*cos(φ)
    u=(-c1,c3)
    v=(-c2,-c4,c5)
    pts = []
    for t in range(0, 2π, length=360)
        pt = c + [u[1]*cos(t) + v[1]*sin(t),u[2]*cos(t)+v[2]*sin(t),v[3]*sin(t)]
        append!(pts,pt)
    end
    return pts
end
end
