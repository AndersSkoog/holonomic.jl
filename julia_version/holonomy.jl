module Holonomy
using LinearAlgebra
using .Primitives:S²,toVec,vec3,vec2,toArray
using .Orientation:SO3
using .Transform:MobiusRotCoef, MobiusTrans
using .Projection:stereographic, inv_stereographic

    function toS²(v::Array{Float})
        x,y,z = v[0],v[1],v[2]
        r = sqrt(x^2+y^2+z^2)
        θ = atan2(y,x)
        φ = acos(z/r)
        return S²(θ,φ)
    end

    struct ArcS²
        from::S²
        to::S²
    end

    function ArcRot(v::ArcS²)
        v1,v1 = toVec(v.from)),toVec(v.to))
        ang = acos(clamp(dot(v1, v2), -1, 1))
        ax  = normalize(cross(v1, v2))
        dir = toS²(ax)
        return SO3(dir,ang)
    end

    function RollTrans(spherecurve::Array{S²}, i::Int, contact::Array{Float,2}, R::SO3)
        r=1.0
        arc = ArcS²(curve[i],curve[i+1])
        ang = acos(clamp(dot(toArray(arc.from), toArray(arc.to)), -1, 1))
        R_inc = rotation(arc)
        R_new = R * R_inc
        n = [R_new[1,3], R_new[2,3], R_new[3,3]]
        d = [0.0,0.0,-1.0]
        move_dir = cross(n, d)
        mdn = norm(move_dir)

        if mdn < 1e-12
            disp = [0.0, 0.0, 0.0]
        else
            move_dir = normalize(move_dir)
            disp = [move_dir[0] * r * ang,
                    move_dir[1] * r * ang,
                    move_dir[2] * r * ang]
        end

        contact_new = contact .+ disp

        return contact_new, R_new
    end

    function inital_R(spherecurve::Array{S²}):
        arc = ArcS²(S²(0,pi/2),curve[0])
        return rotation(arc)
    end

    function SphereCurveRollData(spherecurve::Array{S²})
        prev_contact = vec3(0.0, 0.0, 0.0)
        prev_R = initial_R(curve)
        curveC = vec3[] #planar contact points of rolling sphere
        curveO = SO3[]  #moving frame
        curveT = vec3[] #torsion curve  

        for i in 1:length(curve)-2
            contact, R = rolltranslation(spherecurve,i,prev_contact,prev_R)
            tor_z = R.mtx * [0.0, 0.0, 1.0]
            tor_pt = vec3(contact.x,contact.y,tor_z[3])
            push!(curveC, contact)
            push!(curveO, R)
            push!(curveT, tor_pt)
            prev_contact = contact
            prev_R = R
        end

        return curveC, curveO, curveT
    end

    function cap_circle(sphere_point::Array{Float})
        s2 = toS²(sphere_point)
        n = [0.0,0.0,1.0]
        s = sphere_point
        c = (4/5)*s
        r = sin(pi/5)
        theta=s2.θ
        phi=s2.φ
        c1=r*sin(theta)
        c2=r*sin(phi)*cos(theta)
        c3=r*cos(theta)
        c4=r*sin(phi)*sin(theta)
        c5=r*cos(phi)
        u=(-c1,c3)
        v=(-c2,-c4,c5)
        pts = []
        for t in range(0, 2π, length=360)
          pt = c + [u[0]*cos(t) + v[0]*sin(t),u[1]*cos(t)+v[1]*sin(t),v[2]*sin(t)]
          pts.append!(pt)
        end
        return pts
    end



    function HolonomicView3(sphere_curve::Array{S²},seli::Int)
      fov = sin(pi/5)
      gam = acos(4/5)
      cap_scale = tan(gam/2)
      npole = [0.0,0.0,1.0]
      curveC, curveO, curveT = SphereCurveRollData(sphere_curve)
      sel=Complex(curveC[seli].x,curveC[seli].y)
      a,b,c,d = MobiusRotCoef(sel)
      sphere_point = inv_stereographic_proj(sel)
      persp_point = 2 * sphere_point
      viewpts = []
      transpts = []
      circ = cap_circle(sphere_point)
      p(θ) = cos(γ) * c + sin(γ) * (cos θ * u + sin θ * v)
      for i in 1:length(curve)-2
        p = curveC[i]
        x,y = p[0],p[1]
        scaled_z = Complex(x,y) * cap_scale
        trans_z=MobiusTrans(pz,(a,b,c,d))
        transpts.append!(trans_z)
        vpt = inv_stereographic(trans_z)
        viewpts.append!(vpt)
      end
      return (vp=viewpts,pp=persp_point,sp=sphere_point,tp=transpts)
    end











end







