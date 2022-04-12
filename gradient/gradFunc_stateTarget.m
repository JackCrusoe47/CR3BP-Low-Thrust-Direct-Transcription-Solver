function grad = gradFunc_stateTarget(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;
stateTarget = problem.stateTarget;
weightTarget = problem.weightTarget;

% - gradient of objective
grad = zeros(size(z));
dState = ( z(nState*nSegment*(nOdd-1)+1:nState*nSegment*(nOdd-1)+7) - stateTarget ).*weightTarget;
grad(nState*nSegment*(nOdd-1)+1:nState*nSegment*(nOdd-1)+7) = 2.*dState;

end