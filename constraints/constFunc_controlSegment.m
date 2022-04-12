function Fcontrol = constFunc_controlSegment(zSegment,transcribe,problem)

% - extract parameters
% nState = problem.nState;
nControl = problem.nControl;
% nOdd = transcribe.param.nOdd;
% nEven = transcribe.param.nEven;

% - extract state, control and tof
% xNodes = reshape( zSegment( 1:nState*nOdd ), nState, nOdd );
uSegment = zSegment( end-nControl:end-1 );
% tof = zSegment( end );

% constraint: thrust unit vector uT*uT + uW*uW + uN*nN == 1
Fcontrol = uSegment(2)*uSegment(2) + uSegment(3)*uSegment(3) + uSegment(4)*uSegment(4) - 1;

end