function SDcontrol = strFunc_controlSegment(transcribe,problem)

% SDsegment = [ nState*nEven  x  nState*nOdd + nControl + 1 tof ]

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
% nEven = transcribe.param.nEven;
% nSegment = transcribe.param.nSegment;
timeFixed = problem.flag.timeFixed;

% - extract state, control and tof
% xNodes = reshape( zSegment( 1:nState*nOdd ), nState, nOdd );
% uSegment = zSegment( end-nControl:end-1 );
% tof = zSegment( end );

% initializing
if ~timeFixed
    SDcontrol = zeros( 1 , nState*nOdd + nControl + 1 );
else
    SDcontrol = zeros( 1 , nState*nOdd + nControl );
end
SDcontrol(1,nState*nOdd+2) = 1;
SDcontrol(1,nState*nOdd+3) = 1;
SDcontrol(1,nState*nOdd+4) = 1;

end