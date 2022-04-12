function obj = objFunc_massTarget(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;
massTarget = problem.massTarget;

% - maximize final mass
dMass = ( z(nState*nSegment*(nOdd-1)+7) - massTarget );
obj = dMass ^ 2;

end