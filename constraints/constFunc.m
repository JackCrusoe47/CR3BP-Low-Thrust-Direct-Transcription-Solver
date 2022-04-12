function F = constFunc(z,auxdata)
% -- optimzation constraint

% - extract structs
problem = auxdata.problem;
transcribe = auxdata.transcribe;

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
nSegment = transcribe.param.nSegment;

% - vectors
Fdefect = zeros(nSegment*(nOdd-1)*nState,1);
Fcontrol = zeros(nSegment,1);

% segment initial state and control index
idx0_xNodes = 1;
idx0_uNodes = (nSegment*(nOdd-1)+1)*nState + 1;

% - loop along segments
for i = 1:nSegment
    
    % segment final state and control index
   idxf_xNodes = idx0_xNodes + nOdd*nState - 1;
   idxf_uNodes = idx0_uNodes + nControl - 1;
   
   % segment objective variables
   zSegment = [z(idx0_xNodes:idxf_xNodes);
       z(idx0_uNodes:idxf_uNodes);
       z(end)];
    
   % defect constraint for segment
   Fdefect( (i-1)*(nOdd-1)*nState+1 : (i-1)*(nOdd-1)*nState+(nOdd-1)*nState ) = ...
       constFunc_defectSegment(zSegment,transcribe,problem);
   % control constraint for segment
   Fcontrol(i) = constFunc_controlSegment(zSegment,transcribe,problem);
   
   % update initial indexes
   idx0_xNodes = idxf_xNodes-nState+1;
   idx0_uNodes = idxf_uNodes+1;
   
end

% - path bound constraint
Fpathbnds = constFunc_pathBnds(z,transcribe,problem);

% - full constraint vector
F = [Fdefect; Fcontrol; Fpathbnds];

end