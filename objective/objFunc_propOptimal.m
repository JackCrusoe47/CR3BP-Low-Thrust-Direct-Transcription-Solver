function obj = objFunc_propOptimal(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;

% - maximize final mass
obj = -z(nState*nSegment*(nOdd-1)+7);

end