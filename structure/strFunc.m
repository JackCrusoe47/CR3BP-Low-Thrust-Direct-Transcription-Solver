function SD = strFunc(auxdata)
% -- optimzation jacobian structure

problem = auxdata.problem;
transcribe = auxdata.transcribe;

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
nSegment = transcribe.param.nSegment;
timeFixed = problem.flag.timeFixed;

% - initializing SD matrices

if ~timeFixed
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl + 1;
else
    lenz = (nSegment*(nOdd-1)+1)*nState + nSegment*nControl;
end

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
    
    % - Segment defect
    
    I = (i-1)*(nOdd-1)*nState+1 : (i-1)*(nOdd-1)*nState+(nOdd-1)*nState;
    
    % defect jacobian for segment
    SDdefect_segment = strFunc_defectSegment(transcribe,problem);

    % derivative of defect wrt to segment nodes
    J = idx0_xNodes : idxf_xNodes;
    idxf_defect = idx0_defect + length(I)*length(J) - 1;
    i_defect(idx0_defect:idxf_defect) = reshape( I'.*ones(1,length(J)),length(I)*length(J),1);
    j_defect(idx0_defect:idxf_defect) = reshape( J.*ones(length(I),1),length(I)*length(J),1);    
    v_defect(idx0_defect:idxf_defect) = reshape( SDdefect_segment(:,1:nOdd*nState),length(I)*length(J),1);

    % derivative of defect wrt to segment control
    J = idx0_uNodes : idxf_uNodes;
    idx0_defect = idxf_defect + 1;
    idxf_defect = idx0_defect + length(I)*length(J) - 1;
    i_defect(idx0_defect:idxf_defect) = reshape( I'.*ones(1,length(J)),length(I)*length(J),1);
    j_defect(idx0_defect:idxf_defect) = reshape( J.*ones(length(I),1),length(I)*length(J),1);    
    if ~timeFixed
        v_defect(idx0_defect:idxf_defect) = reshape( SDdefect_segment(:,nOdd*nState+1:end-1),length(I)*length(J),1);
    else
        v_defect(idx0_defect:idxf_defect) = reshape( SDdefect_segment(:,nOdd*nState+1:end),length(I)*length(J),1);
    end
    
    % derivative of defect wrt to time of flight
    if ~timeFixed
        J = lenz;
        idx0_defect = idxf_defect + 1;
        idxf_defect = idx0_defect + length(I)*length(J) - 1;
        i_defect(idx0_defect:idxf_defect) = reshape( I'.*ones(1,1),length(I),1);
        j_defect(idx0_defect:idxf_defect) = reshape( J.*ones(length(I),1),length(I),1);
        v_defect(idx0_defect:idxf_defect) = reshape( SDdefect_segment(:,end),length(I)*length(J),1);
    end
    
    % - Segment control
    
    % control jacobian for segment
    SDcontrol_segment = strFunc_controlSegment(transcribe,problem);
    % derivative of control constraint wrt control
    idxf_control = idx0_control + nControl - 1;
    i_control(idx0_control:idxf_control) = i.*ones(nControl,1);
    j_control(idx0_control:idxf_control) = idx0_uNodes : idxf_uNodes;
    if ~timeFixed
        v_control(idx0_control:idxf_control) = SDcontrol_segment(:,nOdd*nState+1:end-1)';
    else
        v_control(idx0_control:idxf_control) = SDcontrol_segment(:,nOdd*nState+1:end)';
    end
    
    % update initial indexes
    idx0_xNodes = idxf_xNodes-nState+1;
    idx0_uNodes = idxf_uNodes+1;
    idx0_defect = idxf_defect+1;
    idx0_control = idxf_control+1;
end

SDdefect = sparse( i_defect, j_defect, v_defect, nSegment*(nOdd-1)*nState, lenz );
SDcontrol = sparse( i_control, j_control, v_control, nSegment, lenz );

% - Path bounds
SDpathbnds = sparse( strFunc_pathBnds(transcribe,problem) );

% - Full sparse jacobian matrix
SD = [SDdefect;SDcontrol;SDpathbnds];

end