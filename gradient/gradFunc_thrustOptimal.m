function grad = gradFunc_thrustOptimal(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;

% - sum of square of thrust
idx0 = (nSegment*(nOdd-1)+1)*nState + 1;
idxf = (idx0 - 1) + (nSegment-1)*nControl + 1;
idx = idx0:nControl:idxf;

% - gradient of objective
grad = zeros(size(z));
grad(idx) = 2.*z(idx);

end