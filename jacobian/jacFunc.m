function DF = jacFunc(z,auxdata)
% -- optimzation jacobian

problem = auxdata.problem;
transcribe = auxdata.transcribe;

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
nSegment = transcribe.param.nSegment;
timeFixed = problem.flag.timeFixed;

lenz = length(z);

% - initializing DF matrices

% defect matrix vectors and indexes
if ~timeFixed
    n_defect = nSegment*(nOdd-1)*nState*(nOdd*nState+nControl+1);
else
    n_defect = nSegment*(nOdd-1)*nState*(nOdd*nState+nControl);
end
i_defect = zeros(n_defect,1);
j_defect = zeros(n_defect,1);
v_defect = zeros(n_defect,1);
idx0_defect = 1;

% control matrix vectors and indexes
n_control = nSegment*nControl;
i_control = zeros(n_control,1);
j_control = zeros(n_control,1);
v_control = zeros(n_control,1);
idx0_control = 1;

% segment initial state and control index
idx0_xNodes = 1;
idx0_uNodes = (nSegment*(nOdd-1)+1)*nState + 1;

% - loop along segments
for i = 1:nSegment
    
    % segment final state and control index
    idxf_xNodes = idx0_xNodes + nOdd*nState - 1;
    idxf_uNodes = idx0_uNodes + nControl - 1;
    
    % segment objective variables
    
    if ~timeFixed
        zSegment = [z(idx0_xNodes:idxf_xNodes);
        z(idx0_uNodes:idxf_uNodes);
        z(end)];
    else
        zSegment = [z(idx0_xNodes:idxf_xNodes);
        z(idx0_uNodes:idxf_uNodes)];
    end
    
    % - Segment defect
    
    I = (i-1)*(nOdd-1)*nState+1 : (i-1)*(nOdd-1)*nState+(nOdd-1)*nState;
    
    % defect jacobian for segment
    DFdefect_segment = jacFunc_defectSegment(zSegment,transcribe,problem);

    % derivative of defect wrt to segment nodes
    J = idx0_xNodes : idxf_xNodes;
    idxf_defect = idx0_defect + length(I)*length(J) - 1;
    i_defect(idx0_defect:idxf_defect) = reshape( I'.*ones(1,length(J)),length(I)*length(J),1);
    j_defect(idx0_defect:idxf_defect) = reshape( J.*ones(length(I),1),length(I)*length(J),1);    
    v_defect(idx0_defect:idxf_defect) = reshape( DFdefect_segment(:,1:nOdd*nState),length(I)*length(J),1);

    % derivative of defect wrt to segment control
    J = idx0_uNodes : idxf_uNodes;
    idx0_defect = idxf_defect + 1;
    idxf_defect = idx0_defect + length(I)*length(J) - 1;
    i_defect(idx0_defect:idxf_defect) = reshape( I'.*ones(1,length(J)),length(I)*length(J),1);
    j_defect(idx0_defect:idxf_defect) = reshape( J.*ones(length(I),1),length(I)*length(J),1); 
    if ~timeFixed
        v_defect(idx0_defect:idxf_defect) = reshape( DFdefect_segment(:,nOdd*nState+1:end-1),length(I)*length(J),1);
    else
        v_defect(idx0_defect:idxf_defect) = reshape( DFdefect_segment(:,nOdd*nState+1:end),length(I)*length(J),1);
    end
    % derivative of defect wrt to time of flight
    if ~timeFixed
        J = lenz;
        idx0_defect = idxf_defect + 1;
        idxf_defect = idx0_defect + length(I)*length(J) - 1;
        i_defect(idx0_defect:idxf_defect) = reshape( I'.*ones(1,1),length(I),1);
        j_defect(idx0_defect:idxf_defect) = reshape( J.*ones(length(I),1),length(I),1);
        v_defect(idx0_defect:idxf_defect) = reshape( DFdefect_segment(:,end),length(I)*length(J),1);
    end
    
    % - Segment control
    
    % control jacobian for segment
    DFcontrol_segment = jacFunc_controlSegment(zSegment,transcribe,problem);
    % derivative of control constraint wrt control
    idxf_control = idx0_control + nControl - 1;
    i_control(idx0_control:idxf_control) = i.*ones(nControl,1);
    j_control(idx0_control:idxf_control) = idx0_uNodes : idxf_uNodes;
    if ~timeFixed
        v_control(idx0_control:idxf_control) = DFcontrol_segment(:,nOdd*nState+1:end-1)';
    else
        v_control(idx0_control:idxf_control) = DFcontrol_segment(:,nOdd*nState+1:end)';
    end
    
    % update initial indexes
    idx0_xNodes = idxf_xNodes-nState+1;
    idx0_uNodes = idxf_uNodes+1;
    idx0_defect = idxf_defect+1;
    idx0_control = idxf_control+1;
end

DFdefect = sparse( i_defect, j_defect, v_defect, nSegment*(nOdd-1)*nState, lenz );
DFcontrol = sparse( i_control, j_control, v_control, nSegment, lenz );

% - Path bounds
DFpathbnds = sparse( jacFunc_pathBnds(z,transcribe,problem) );

% - Full sparse jacobian matrix
DF = [DFdefect;DFcontrol;DFpathbnds];

end