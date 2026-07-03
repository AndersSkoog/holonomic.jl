module H

    struct UnitQuaternion <: FieldVector{4,Float64}
        w::Float64
        i::Float64
        j::Float64
        k::Float64
    end


    function UnitQuaternion(axis::S2, ang::Float64)
        s = sin(ang / 2)
        w = cos(ang / 2)
        i = c.x * s
        j = c.y * s
        k = c.z * s
        n = sqrt(w*w + i*i + j*j + k*k)
        UnitQuaternion(w/n, i/n, j/n, k/n)
    end
end

