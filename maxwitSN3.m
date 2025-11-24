function alpha = maxwitSN3(r)

%-------------------------------------------------------------------------%
%This function computes the upper bounds on the Schmidt number witness for
%nc = 3 copies.
%
%Input: 
% - r: Schmidt number
%Output:
% - alpha: upper bound on the Schmidt number witness Wn.
%-------------------------------------------------------------------------%

d = 2; %Local dimension
nc = 3; %Number of copies
id = eye(d); %Local identity

%Computational basis
for l = 0 : d-1
    comp{l+1} = id(:,l+1);
end

%Computational basis of A1A2A3 (or B1B2B3)
id123 = Tensor(id,nc);

for l = 0 : d^nc-1
    basis123{l+1} = id123(:,l+1);
end

%Computational basis on A' (or B')
idr = eye(r); %Schmidt rank identity

for l = 0 : r-1
    basisprime{l+1} = idr(:,l+1);
end

%Enlarging state
phi = 0;
for l = 0 : r-1
    phi = phi + 1/sqrt(r)*Tensor(idr(:,l+1),2);
end

idnc = Tensor(id,2*nc); %Identity of A1A2A3 B1B2B3

%Pauli strings witness
G = PauliStrings(nc);
W = 0;
for k = 1 : 3*2+1
    W = W + Tensor(G{nc,k},transpose(G{nc,k}));
end

%Enlarged witness
WL = Tensor(W,phi*phi');

%Permutation of subsystems
[~,perm] = sort([1 2 3 5 6 7 4 8]); %A1A2A3 B1B2B3 A'B' -> A1A2A3A' B1B2B3B'
WL = PermuteSystems(WL,perm,[d d d d d d r r]);

%% SDP

alpha1 = 1;
alpha2 = 0;

%Initial B state
sigmaBBb = RandomDensityMatrix(d^nc*r);

while abs(alpha1-alpha2) > 10e-8

    alpha1 = alpha2;

    %Variables
    sigmaAAa = sdpvar(d^nc*r);

    %Constraints
    C = [sigmaAAa >= 0, trace(sigmaAAa) == 1];

    %Enlarged state
    psi = Tensor(sigmaAAa,sigmaBBb);

    %SolveSDP
    disp('Options')
    ops=sdpsettings('solver','mosek', 'cachesolvers', 1);
    diagnostic=solvesdp(C,-real(trace(WL*psi)),ops)

    sigmaAAa = double(sigmaAAa);


    %Variables
    sigmaBBb = sdpvar(d^nc*r);

    %Constraints
    C = [sigmaBBb >= 0, trace(sigmaBBb) == 1, PartialTrace(sigmaBBb,1,[d^nc,r]) == idr/r];

    %Enlarged separable state
    psi = Tensor(sigmaAAa,sigmaBBb);

    %SolveSDP
    disp('Options')
    ops=sdpsettings('solver','mosek', 'cachesolvers', 1);
    diagnostic=solvesdp(C,-real(trace(WL*psi)),ops)

    sigmaBBb = double(sigmaBBb);

    alpha2 = r^2*double(trace(WL*psi));

end

alpha = alpha2;

end