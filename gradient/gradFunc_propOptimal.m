function grad = gradFunc_propOptimal(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;

% - gradient of objective
grad = zeros(size(z));
grad(nState*nSegment*(nOdd-1)+7) = -1;

end