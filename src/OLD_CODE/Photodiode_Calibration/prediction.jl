function f(x,w)
    # State Propagation only involves a slight rotation 
    q = x[1:4];     # Quaternion portion
    b = x[5:7];     # Bias portion

    γ = w-b;        # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    theta = (nγ*dt);  
    r = γ/nγ;  # Make unit

    x[1:4] = qmult(q, [r*sin(theta/2); cos(theta/2)]); 
    return x
end

function prediction(xk, w, dt, numDiodes)
    # Predicts next state using current state and angular velocity
    #
    # Arguments:
    #   - xk: Current state [q β c α ϵ]                  | [7 + 3i,]
    #   - w: Current angular velocity                    | [3,]
    #   - dt: Time step                                  | Scalar
    #   - numDiodes: Number of photodiodes being used    | Scalar
    #
    # Returns:
    #   - xn: Predicted next state [q β c α ϵ]                     | [7 + 3i,]
    #   - A: State Jacobian with respect to state    
    #           (note that quaternions are replaced with 3 param)  | [(6 + 3i) x (6 + 3i)]

    q = xk[1:4]; # Quaternion portion
    b = xk[5:7]; # Bias portion

    γ = w-b;     # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    theta = (nγ*dt);  
    r = γ/nγ;  # Make unit

    qp = qmult(q, [r*sin(theta/2); cos(theta/2)]); 
    
    skew = -1*[0 -γ[3] γ[2]; γ[3] 0 -γ[1]; -γ[2] γ[1] 0]; # Hat

    R = (I(3) + (skew/nγ)*sin(nγ*dt) + ((skew/nγ)^2)*(1 - cos(nγ*dt)));     # Rodrigues (for matrix exponential?)

    A = [R -dt*I(3); zeros(3,3) I(3)]; # Jacobian of f(x)

    c = xk[8:(7+numDiodes)]
    α = xk[(8+numDiodes):(7 + 2*numDiodes)]
    ϵ = xk[(8 + 2*numDiodes):end]

    A = [A                           zeros(6, 3 * numDiodes); 
         zeros(3 * numDiodes, 6)     I(3*numDiodes)];

    xn = [qp; b; c; α; ϵ]; # x at next step

    return xn, A
end