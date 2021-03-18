function triad(rN1,rN2,rB1,rB2)
    # Method for estimating the rotation matrix between two reference frames. 
    #   - Takes in a pair of vectors in the newtonian (inertial)
    #     frame and their corresponding pairs in the body frame.
    #   - Generates a DCM and then converts it into a quaternion (scalar
    #       first)

    tN1 = rN1;
    tN2 = cross(rN1,rN2)/norm(cross(rN1,rN2));
    tN3 = cross(tN1,tN2);

    nT = [tN1[:] tN2[:] tN3[:]];


    tB1 = rB1;
    tB2 = cross(rB1,rB2)/norm(cross(rB1,rB2));
    tB3 = cross(tB1,tB2);

    bT = [tB1[:] tB2[:] tB3[:]];

    # DCM
    R = nT*(bT');

    # QUATERNION
    q = q_from_DCM(R);

    return q, R
end
