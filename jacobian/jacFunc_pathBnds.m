function DFpathbnds = jacFunc_pathBnds(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;

% - path end point state bounds
DFpathbnds = zeros(2*nState,length(z));
DFpathbnds(1:nState,1:nState) = eye(nState);
DFpathbnds(nState+1:2*nState,nState*(nOdd-1)*nSegment+1:nState*(nOdd-1)*nSegment+nState) = eye(nState);

end