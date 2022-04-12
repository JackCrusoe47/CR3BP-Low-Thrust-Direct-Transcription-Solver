function SDpathbnds = strFunc_pathBnds(transcribe,problem)

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;
timeFixed = problem.flag.timeFixed;

% - initializing SD matrices

if ~timeFixed
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl + 1;
else
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl;
end

% - path end point state bounds
SDpathbnds = zeros(2*nState,lenz);
SDpathbnds(1:nState,1:nState) = eye(nState);
SDpathbnds(nState+1:2*nState,nState*(nOdd-1)*nSegment+1:nState*(nOdd-1)*nSegment+nState) = eye(nState);

end