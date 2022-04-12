function DFcontrol = jacFunc_controlSegment(zSegment,transcribe,problem)

% DFsegment = [ nState*nEven  x  nState*nOdd + nControl + 1 tof ]

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
% nEven = transcribe.param.nEven;
% nSegment = transcribe.param.nSegment;
timeFixed = problem.flag.timeFixed;

% - extract state, control and tof
% xNodes = reshape( zSegment( 1:nState*nOdd ), nState, nOdd );
uSegment = zSegment( end-nControl:end-1 );
% tof = zSegment( end );

% initializing
if ~timeFixed
    DFcontrol = zeros( 1 , nState*nOdd + nControl + 1 );
else
    DFcontrol = zeros( 1 , nState*nOdd + nControl );
end
DFcontrol(1,nState*nOdd+2) = 2*uSegment(2);
DFcontrol(1,nState*nOdd+3) = 2*uSegment(3);
DFcontrol(1,nState*nOdd+4) = 2*uSegment(4);

end