module Orientation

    using .Primitives:S²,vec2,vec3,toVec,toArray
    export SO3, UnitQuaternion, SU2, MobiusRotCoef

    # -------------------------
    # SO(3)
    # -------------------------
    struct SO3
        mtx::Matrix{Float}

        function SO3(dir::S²,ang::Float)
            x, y, z = vec3(dir)
            c = cos(ang)
            s = sin(ang)
            C = 1 - c
            R = [
                c + x*x*C   x*y*C - z*s   x*z*C + y*s
                y*x*C + z*s  c + y*y*C     y*z*C - x*s
                z*x*C - y*s  z*y*C + x*s   c + z*z*C
            ]
            new(R)
        end
    end

    # -------------------------
    # Unit Quaternion
    # -------------------------
    struct UnitQuaternion
        w::Float
        i::Float
        j::Float
        k::Float

        function UnitQuaternion(dir::S², ang::Float)
            s = sin(ang / 2)
            w = cos(ang / 2)
            c = toVec(dir)

            i = c.x * s
            j = c.y * s
            k = c.z * s

            n = sqrt(w*w + i*i + j*j + k*k)

            new(w/n, i/n, j/n, k/n)
        end
    end


    # -------------------------
    # SU(2) / S³ representation
    # -------------------------
    struct S³
        U::Matrix{Complex}

        function S³(dir::S², ang::Float)
            c = toVec(dir)
            σ1 = [0+0im 1+0im; 1+0im 0+0im]
            σ2 = [0+0im 0-1im; 0+1im 0+0im]
            σ3 = [1+0im 0+0im; 0+0im -1+0im]
            I  = [1+0im 0+0im; 0+0im 1+0im]
            dotsum = c.x*σ1 + c.y*σ2 + c.z*σ3
            U = cos(ang/2)*I - im*sin(ang/2)*dotsum
            new(U)
        end

    end

end
