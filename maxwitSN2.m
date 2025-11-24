function alpha = maxwitSN2(r)

%-------------------------------------------------------------------------%
%This function computes the upper bounds on the Schmidt number witness for
%nc = 2 copies.
%
%Input: 
% - r: Schmidt number
%Output:
% - alpha: upper bound on the Schmidt number witness Wn.
%-------------------------------------------------------------------------%


d = 2; %Local dimension
nc = 2; %Number of copies
id = eye(d); %Local identity

%Computational basis
for l = 0 : d-1
    comp{l+1} = id(:,l+1);
end

%Computational basis of A1A2 (or B1B2)
id12 = Tensor(id,nc);

for l = 0 : d^nc-1
    basis12{l+1} = id12(:,l+1);
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

idnc = Tensor(id,2*nc); %Identity of A1A2 B1B2

%Pauli strings witness
G = PauliStrings(nc);
W = 0;
for k = 1 : 2*nc+1
    W = W + Tensor(G{nc,k},transpose(G{nc,k}));
end

%Enlarged witness
WL = Tensor(W,phi*phi');

%Permutation of subsystems
[~,perm] = sort([1 2 4 5 3 6]); %A1A2 B1B2 A'B' -> A1A2A' B1B2B'
WL = PermuteSystems(WL,perm,[d d d d r r]);

%% SDP

%Variables
sigmaAAa = sdpvar(d^nc*r);

%Constraints
C = [sigmaAAa >= 0, trace(sigmaAAa) == 1];

%B state
phiBBb = 0;
for k = 1 : r
    phiBBb = phiBBb + 1/sqrt(r)*Tensor(id12(:,k),idr(:,k));
end
sigmaBBb = phiBBb*phiBBb';

%Enlarged state
psi = Tensor(sigmaAAa,sigmaBBb);

%SolveSDP
disp('Options')
ops=sdpsettings('solver','mosek', 'cachesolvers', 1);
diagnostic=solvesdp(C,-real(trace(WL*psi)),ops)

alpha = r^2*double(trace(WL*psi));


end