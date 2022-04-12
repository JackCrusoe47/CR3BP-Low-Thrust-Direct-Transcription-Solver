function SDsegment = strFunc_defectSegment(transcribe,problem)

% SDsegment = [ nState*nEven  x  nState*nOdd + nControl + 1 tof ]

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
nEven = transcribe.param.nEven;
% nSegment = transcribe.param.nSegment;
timeFixed = problem.flag.timeFixed;

% the jacobian is almost fully populated
if ~timeFixed
    SDsegment = ones( nState*nEven , nState*nOdd + nControl + 1 );
else
    SDsegment = ones( nState*nEven , nState*nOdd + nControl );
end

end