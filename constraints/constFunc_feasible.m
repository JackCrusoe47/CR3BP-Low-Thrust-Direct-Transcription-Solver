function F = constFunc_feasible(z,auxdata)
% -- optimzation constraint

% - extract structs
problem = auxdata.problem;
transcribe = auxdata.transcribe;

nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
nEven = transcribe.param.nEven;
nSegment = transcribe.param.nSegment;
timeFixed = problem.flag.timeFixed;
maxThrust = problem.param.maxThrust;

flag0 = problem.flag.path0;
flagf = problem.flag.pathf;

path0 = problem.bnds.path0.ub;
pathf = problem.bnds.pathf.ub;

if ~timeFixed
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl + 1;
else
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl;
end

% - compute transcription constraints
F = constFunc(z(1:lenz),auxdata);

% - seperate transcription constraints
Fdefect = F(1:(nSegment*nEven)*nState);
Fcontrol = F((nSegment*nEven)*nState+1:(nSegment*nEven)*nState+nSegment);
Fpathbnds = F((nSegment*nEven)*nState+nSegment+1:end);

% remove unecessary path boundary constraints
Fpathbnds = Fpathbnds - [path0;pathf];
Fpathbnds = Fpathbnds( [flag0;flagf] == 1 );

idx0_uNodes = (nSegment*(nOdd-1)+1)*nState + 1;
idxf_uNodes = (nSegment*(nOdd-1)+1)*nState + (nSegment-1)*nControl + 1;

% extract thrust and slack variables
zThrust = z( idx0_uNodes : nControl : idxf_uNodes );
zSlack = z(lenz+1:end);

% add thrust inequality constraints
Fthrust = zThrust - maxThrust.*sin(zSlack).^2;

% modified constraint vector
F = [Fdefect; Fcontrol; Fpathbnds; Fthrust];

end