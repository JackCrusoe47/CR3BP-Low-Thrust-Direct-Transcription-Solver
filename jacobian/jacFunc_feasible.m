function DF = jacFunc_feasible(z,auxdata)

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

if ~timeFixed
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl + 1;
else
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl;
end

% -- compute regular jacobian
DF = jacFunc(z(1:lenz),auxdata);

% - seperate transcription constraints
DFdefect = DF(1:(nSegment*nEven)*nState,:);
DFcontrol = DF((nSegment*nEven)*nState+1:(nSegment*nEven)*nState+nSegment,:);
DFpathbnds = DF((nSegment*nEven)*nState+nSegment+1:end,:);

% - remove unecessary jacobian rows
DFpathbnds = DFpathbnds( [flag0;flagf] == 1 , : );

DF = [DFdefect;DFcontrol;DFpathbnds];

% - add zero right rectangular matrix
DFright = sparse([], [], [], size(DF,1), nSegment);

% - bottom rectangular matrix
DFbottom1 = sparse([], [], [], nSegment, (nSegment*(nOdd-1)+1)*nState);
DFbottom2 = sparse(1:nSegment, 1:nControl:nControl*nSegment, ones(1,nSegment), nSegment, nSegment*nControl);
if ~timeFixed
    DFbottom3 =  sparse([], [], [], nSegment, 1 );
    DFbottom = [ DFbottom1, DFbottom2, DFbottom3 ];
else
    DFbottom = [ DFbottom1, DFbottom2 ];
end
    

% - slack varaibles
zSlack = z(lenz+1:end);

% - add thrust constraint jacobain bottom right matrix
% derivative of u - umax*sin(s)^2  = -2*umax*sin(s)*cos(s)
DFthrust = speye(nSegment) .* (-2*maxThrust.*sin(zSlack).*cos(zSlack));

% - modified jacobian matrix
DF = [ DF, DFright;
    DFbottom, DFthrust ];

end