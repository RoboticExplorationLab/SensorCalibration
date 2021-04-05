function f(x,w)
    # State Propagation only involves a slight rotation 
    q = x[1:4];     # Quaternion portion
    b = x[5:7]; # Bias portion

    γ = w-b;        # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    theta = (nγ*dt);  
    r = γ/nγ;  # Make unit

    x[1:4] = qmult(q, [r*sin(theta/2); cos(theta/2)]); 
    return x
end

function prediction(xk,w,dt, numDiodes)
    q = xk[1:4]; # Quaternion portion
    b = xk[5:7]; # Bias portion

    γ = w-b;        # Adjusted angular velocity (w - biases)
    nγ = norm(γ)

    theta = (nγ*dt);  
    r = γ/nγ;  # Make unit

    qp = qmult(q, [r*sin(theta/2); cos(theta/2)]); 
    
    skew = -1*[0 -γ[3] γ[2]; γ[3] 0 -γ[1]; -γ[2] γ[1] 0]; # Hat

    R = (I(3) + (skew/nγ)*sin(nγ*dt) + ((skew/nγ)^2)*(1 - cos(nγ*dt)));     # Rodrigues 

    # A is jacobian of f(x)
    A = [R -dt*I(3); zeros(3,3) I(3)];
    # xn = [qp; b];

    c = xk[8:(7+numDiodes)]
    α = xk[(8+numDiodes):(7 + 2*numDiodes)]
    ϵ = xk[(8 + 2*numDiodes):end]

    A = [A                           zeros(6, 3 * numDiodes); 
         zeros(3 * numDiodes, 6)     I(3*numDiodes)];

    # _fClosure(x) = f(x, w);
    # A_alt = ForwardDiff.jacobian(_fClosure, xk)
    # println("Size A: ", size(A))
    # println("Size A2: ", size(A_alt))
    # diff = abs.(A - A_alt);
    # println("Diff: ", sum(diff))

    xn = [qp; b; c; α; ϵ]; # x at next step

    return xn, A

end
