function triad(rN1,rN2,rB1,rB2)
    # Method for estimating the rotation matrix between two reference frames 
    #   Relies only on a single pair of vectors in each frame
    # Inputs: 
    #   - rN1, rN2: Pair of vectors in the Newtonian (inertial) frame     | [3,]
    #   - rB1, rB2: Corresponding pair of vectors in body frame           | [3,]    # 
    # Outputs:
    #   - R: A directed cosine matrix (DCM) representing the rotation     | [3 x 3]
    #           between the two frames 
    #   - q: A quaternion (scalar last) representing the rotation         | [4,]
    #           between the two frames


    # This method avoids taking an inverse by generating an orthogonal matrix
    tN1 = rN1;
    tN2 = cross(rN1,rN2)/norm(cross(rN1,rN2));
    tN3 = cross(tN1,tN2)/norm(cross(tN1,tN2));  

    nT = [tN1[:] tN2[:] tN3[:]];


    tB1 = rB1;
    tB2 = cross(rB1,rB2)/norm(cross(rB1,rB2));
    tB3 = cross(tB1,tB2)/norm(cross(tB1,tB2))

    bT = [tB1[:] tB2[:] tB3[:]];

    # DCM
    R = nT*(bT');

    # QUATERNION
    q = q_from_DCM(R);

    return q, R
end
