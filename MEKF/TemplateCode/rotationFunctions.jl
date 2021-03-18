function q_from_DCM(R)
    # Converts a scalar-first quaternion to a DCM matrix (From Kevin)

    T = R[1,1] + R[2,2] + R[3,3];
    if (T > R[1,1]) && (T > R[2,2]) && (T>R[3,3])
        q4 = .5*sqrt(1+T);
        r  = .25/q4;
        q1 = (R[3,2] - R[2,3])*r;
        q2 = (R[1,3] - R[3,1])*r;
        q3 = (R[2,1] - R[1,2])*r;
    elseif R[1,1]>R[2,2] && R[1,1]>R[3,3]
        q1 = .5*sqrt(1-T + 2*R[1,1]);
        r  = .25/q1;
        q4 = (R[3,2] - R[2,3])*r;
        q2 = (R[1,2] + R[2,1])*r;
        q3 = (R[1,3] + R[3,1])*r;
    elseif R[2,2]>R[3,3]
        q2 = .5*sqrt(1-T + 2*R[2,2]);
        r  = .25/q2;
        q4 = (R[1,3] - R[3,1])*r;
        q1 = (R[1,2] + R[2,1])*r;
        q3 = (R[2,3] + R[3,2])*r;
    else
        q3 = .5*sqrt(1-T + 2*R[3,3]);
        r  = .25/q3;
        q4 = (R[2,1] - R[1,2])*r;
        q1 = (R[1,3] + R[3,1])*r;
        q2 = (R[2,3] + R[3,2])*r;
    end
    q = [q1;q2;q3;q4];     # Scalar last ##########################
    # q = [q4;q1;q2;q3];     # Scalar first
    if q4<0
       q = -q;
    end

    return q
end

function dcm_from_q(quat)
    # Takes in a quaternion (scalar first - [q, q⃗]) and returns the DCM  (From Kevin)
    
    quat = quat[:];
    
    # v = quat[2:4];                    #############
    # s = quat[1];
    v = quat[1:3];
    s = quat[4];
    
    DCM = I(3) +2*hat(v)*(s*I(3) + hat(v));

    return DCM;
    
end

function hat(v)
    # Forms a skew symmetric matrix from a vector 

    v = v[:];
    vhat = [  0   -v[3]  v[2]; 
             v[3]   0   -v[1];
            -v[2]  v[1]   0];

    return vhat;
end

function qmult(q, p)
    # Performs quaternion multiplication (Hamilton product) 
    #   (Quaternions are scalar first)

    q = q[:]
    p = p[:]

    # q0, p0 = q[1], p[1]           ##############33
    # q⃗ = q[2:4]
    # p⃗ = p[2:4]

    # r0 = q0 * p0 - dot(q⃗, p⃗)
    # r⃗ = (q0 * p⃗) + (p0 * q⃗) + cross(q⃗, p⃗)

    # r = [r0; r⃗]
    
    q0, p0 = q[4], p[4]
    q⃗ = q[1:3]
    p⃗ = p[1:3]

    r0 = q0 * p0 - dot(q⃗, p⃗)
    r⃗ = (q0 * p⃗) + (p0 * q⃗) + cross(q⃗, p⃗)

    r = [r⃗; r0]





    return r
end

function qconj(q) ##############
    # Returns the quaternion conjugate (scalar first)
    # q0 = q[1]
    # q⃗ = q[2:4]

    # return [q0; -q⃗]

    q0 = q[4]
    q⃗ = q[1:3]
    return [-q⃗; q0]

end
    