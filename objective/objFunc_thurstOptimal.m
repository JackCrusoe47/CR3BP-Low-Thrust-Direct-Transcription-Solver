function obj = objFunc_thurstOptimal(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;

% - sum of square of thrust
idx0 = (nSegment*(nOdd-1)+1)*nState + 1;
idxf = (idx0 - 1) + (nSegment-1)*nControl + 1;
idx = idx0:nControl:idxf;
obj = sum( z(idx).^2 );

end