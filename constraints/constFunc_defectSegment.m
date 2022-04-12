function Fdefect = constFunc_defectSegment(zSegment,transcribe,problem)

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
nEven = transcribe.param.nEven;
timeFixed = problem.flag.timeFixed;

% - extract state, control and tof
xNodes = reshape( zSegment( 1:nState*nOdd ), nState, nOdd );
uSegment = zSegment( end-nControl:end-1 );
if ~timeFixed
    tof = zSegment( end );
else
    tof = transcribe.traj.nodes.t(end);
end

% - trascription matrices and vectors
Ainv = transcribe.param.Ainv;
B = transcribe.param.B;
D = transcribe.param.D;
W = transcribe.param.weightEven;

% - function handles
dyn = problem.func.dynamics;

% - segment length
Deltat = tof/transcribe.param.nSegment;

% - constraint vector
Fdefect = zeros(7*nEven,1);

% matrix of node states and derivatives
X = [ xNodes, Deltat/2*dyn(ones(1,nOdd), xNodes, uSegment*ones(1,nOdd) ) ];

% coefficients of nth order polynomial for states
C = X*Ainv;

% compute defect constraints at even points
for i = 1:nEven
    % ith even defect constraint (weighted for performance)
    Fdefect( (i-1)*nState+1 : (i-1)*nState+nState, 1 ) = ...
        ( C*D(:,i) - Deltat/2*dyn( 1, C*B(:,i), uSegment ) )*W(i);
end

end