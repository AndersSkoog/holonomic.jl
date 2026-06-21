module Primitives

    export S¹, B², S², B³, C², vec2, vec3, vec4, toVec, toArray

    struct S¹
        θ::Float
        function S¹(θ::Float)
            θ = mod(θ, 2π)
            new(θ)
        end
    end

    struct B²
        r::Float
        θ::Float
        function B²(θ::Float, r::Float)
            r = clamp(r, 0.0, 1.0)
            θ = mod(θ, 2π)
            new(θ,r)
        end
    end

    struct S²
        θ::Float
        φ::Float
        function S²(θ::Float, φ::Float)
            new(mod(θ, 2π), mod(φ, 2π))
        end
    end

    struct B³
        r::Float
        θ::Float
        φ::Float
        function B³(θ::Float, φ::Float, r::Float)
            new(clamp(r, 0.0, 1.0),mod(θ, 2π),mod(φ, 2π))
        end
    end

    struct C²
        ζ::Complex
        ξ::Complex    
    end

    struct vec2
        x::Float
        y::Float
    end

    struct vec3
        x::Float
        y::Float
        z::Float
    end

    struct vec4
        x::Float
        y::Float
        z::Float
        w::Float
    end

    function toS²(v::vec3)
        r = sqrt(v.x^2+v.y^2+v.z^2)
        θ = atan2(v.y,v.x)
        φ = acos(v.z / r)
        return S²(θ,φ)
    end

    toVec(v::S¹) = vec2(sin(v.θ),cos(v.θ))
    toVec(v::S²) = vec2(sin(v.θ),cos(v.θ))
    toVec(v::B²) = vec2(v.r * sin(v.θ),v.r * cos(v.θ))
    toVec(v::S²) = vec3(cos(v.φ) * cos(v.θ),cos(v.φ) * sin(v.θ),sin(v.φ))
    toVec(v::B³) = vec3(v.r * (cos(v.φ) * cos(v.θ)),v.r * (cos(v.φ) * sin(v.θ)),v.r * sin(v.φ))
    toVec(v::Complex) = vec2(real(v),imag(v))

    toArray(v::vec3) = [v.x,v.y,v.z])
    toArray(v::vec2) = [v.x,v.y])
    toArray(v::vec4) = [v.x,v.y,v.z,v.w])
    toArray(v::S¹) = toArray(toVec(v))
    toArray(v::B²) = toArray(toVec(v))
    toArray(v::S²) = toArray(toVec(v))
    toArray(v::B³) = toArray(toVec(v))
    toArray(v::Complex) = toArray(toVec(v))
end

