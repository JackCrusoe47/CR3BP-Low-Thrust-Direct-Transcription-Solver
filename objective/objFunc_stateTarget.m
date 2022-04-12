function obj = objFunc_stateTarget(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;
stateTarget = problem.stateTarget;
weightTarget = problem.weightTarget;

% - maximize final mass
dState = ( z(nState*nSegment*(nOdd-1)+1:nState*nSegment*(nOdd-1)+7) - stateTarget ).*weightTarget;
obj = sum( dState .^ 2 );

end