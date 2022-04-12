function grad = gradFunc_massTarget(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;
massTarget = problem.massTarget;

% - gradient of objective
grad = zeros(size(z));
dMass = ( z(nState*nSegment*(nOdd-1)+7) - massTarget );
grad(nState*nSegment*(nOdd-1)+7) = 2*dMass;

end