function paramTranscribe = getParamTranscribe_LGL(nOrder)

% -- compute LGL points

% - sanity checks

% minimum order check (minimum order for Hermit Simpson)
if nOrder < 3
    fprintf(2,"WARNING: Interpolation order < 3. Setting to 3!...\n");
    nOrder = 3;
end

% even order check (accuracy same as nearest odd order)
if mod(nOrder,2) == 0
    fprintf(2,"WARNING: Interpolation order even. Setting nearest odd!...\n");
    nOrder = nOrder - 1;
end

% - get LGL points
[tau,weight] = nodes_LGL(nOrder);

% odd and even points
tauOdd = tau(1:2:nOrder);
tauEven = tau(2:2:nOrder-1);

% Number of odd and event points
nOdd = (nOrder+1)/2;
nEven = (nOrder-1)/2;

% odd and even wieghts
weightOdd = weight(1:2:nOrder);
weightEven = weight(2:2:nOrder-1);

% - save parameters
paramTranscribe.nOrder = nOrder;        % order of direct trascription
paramTranscribe.tau = tau;              % LGL node points
paramTranscribe.weight = weight;        % LGL node wieghts
paramTranscribe.nOdd = nOdd;            % LGL number of odd points
paramTranscribe.nEven = nEven;          % LGL number of even points
paramTranscribe.tauOdd = tauOdd;        % LGL odd points
paramTranscribe.tauEven = tauEven;      % LGL even points
paramTranscribe.weightOdd = weightOdd;  % LGL odd weights
paramTranscribe.weightEven = weightEven;% LGL odd weights


% -- compute A matrix

% - submatrix A1 (tau matrix)
A1 = zeros(nOrder+1,nOdd);
A1(1,:) = 1;
for i = 1:nOrder
   A1(1+i,:) = tauOdd.^i; 
end

% - submatrix A2 (dtau matrix)
A2 = zeros(nOrder+1,nOdd);
for i = 1:nOrder
   A2(1+i,:) = i*tauOdd.^(i-1); 
end

% - matrix A
A = [A1,A2];
Ainv = inv(A);

% -- compute B matrix
B = zeros(nOrder+1,nEven);
B(1,:) = 1;
for i = 1:nOrder
   B(1+i,:) = tauEven.^(i); 
end

% -- compute D matrix
D = zeros(nOrder+1,nEven);
D(1,:) = 0;
for i = 1:nOrder
   D(1+i,:) = i*tauEven.^(i-1); 
end

% -- save matrices
paramTranscribe.A = A;
paramTranscribe.Ainv = Ainv;
paramTranscribe.B = B;
paramTranscribe.D = D;

end