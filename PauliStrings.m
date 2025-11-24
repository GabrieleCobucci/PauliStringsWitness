function G = PauliStrings(nc)

%-------------------------------------------------------------------------%
%This function computes the Pauli strings for a given number of copies.
%
%Input: 
% - d: Local dimension
% - nc: Number of copies;
%Output:
% - G{l,:}: Set of Pauli strings for l = 1 : nc
%-------------------------------------------------------------------------%

d = 2; %Local dimension
id = eye(d);

Z = GenPauli(0,1,d); %Pauli Z
X = GenPauli(1,0,d); %Pauli X
Y = i*X*Z; %Pauli Y

G{1,1} = Z;
G{1,2} = X;
G{1,3} = Y;

for l = 2 : nc
    for k = 1 : 2*l - 1
        G{l,k} = Tensor(G{l-1,k},Z);
    end
    G{l,2*l} = Tensor(Tensor(id,l-1),X);
    G{l,2*l+1} = Tensor(Tensor(id,l-1),Y);
end

end