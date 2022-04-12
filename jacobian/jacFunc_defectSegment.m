function DFsegment = jacFunc_defectSegment(zSegment,transcribe,problem)

% DFsegment = [ nState*nEven  x  nState*nOdd + nControl + 1 tof ]

% - extract parameters
nState = problem.nState;
nControl = problem.nControl;
nOdd = transcribe.param.nOdd;
nEven = transcribe.param.nEven;
nSegment = transcribe.param.nSegment;
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
jac = problem.func.jacobian;

% - segment length
Deltat = tof/transcribe.param.nSegment;

% initializing DF vector
if ~timeFixed
    DFsegment = zeros( nState*nEven , nState*nOdd + nControl + 1 );
else
    DFsegment = zeros( nState*nEven , nState*nOdd + nControl );
end

% matrix of node states and derivatives
X = [ xNodes, Deltat/2*dyn(ones(1,nOdd), xNodes, uSegment*ones(1,nOdd) ) ];

% coefficients of nth order polynomial for states
C = X*Ainv;

% - precompute jacobians
% jacobian of node points
dxGradOdd = zeros(nState,nState+nControl);
for j = 1:nOdd
    % jacobian of dynamics wrt odd state and control
    dxGradOdd(:,:,j) = jac(1,xNodes(:,j),uSegment);
end
% jacobian of defect points
dxGradEven = zeros(nState,nState+nControl);
for i = 1:nEven
    % jacobian of dynamics wrt even state and control
    dxGradEven(:,:,i) = jac(1,C*B(:,i),uSegment);
end

% compute defect constraints at even points
for i = 1:nEven
    
    % -- Compute derivative wrt odd nodes
    
    % loop around odd nodes
    for j = 1:nOdd
        
        for k = 1:nState
            
            % derivative of X with kth state of jth node
            dX = zeros(nState,2*nOdd);
            dX(k,j) = 1;
            dX(:,nOdd+j) = Deltat/2*dxGradOdd(:,k,j);
            
            % derivative of defect i with node j
            DFsegment( (i-1)*nState+1 : (i-1)*nState+nState, (j-1)*nState+k  ) = ...
                (dX*Ainv*D(:,i) - ...
                Deltat/2*dxGradEven(:,1:nState,i)*dX*Ainv*B(:,i))*W(i);
            
        end
    end
    
    % -- Compute derivative wrt control
    
    % loop around each control
    for k = 1:nControl
        % derivative of X with kth control
        dX = zeros(nState,2*nOdd);
        for j = 1:nOdd
            dX(:,nOdd+j) = Deltat/2*dxGradOdd(:,nState+k,j);
        end
        
        % derivative of defect i with node j
            DFsegment( (i-1)*nState+1 : (i-1)*nState+nState, nOdd*nState+k  ) = ...
                (dX*Ainv*D(:,i) - ...
                Deltat/2*dxGradEven(:,1:nState,i)*dX*Ainv*B(:,i) - ...
                Deltat/2*dxGradEven(:,nState+k,i))*W(i);
        
    end
    
    % -- Compute derivative wrt time of flight
    
    if ~timeFixed
    
        % derivative of X with tof
        dX = zeros(nState,2*nOdd);
        for j = 1:nOdd
            dX(:,nOdd+j) = 1/2*1/nSegment*dyn(1,xNodes(:,j),uSegment);
        end
        
        % derivative of defect i with tof
        DFsegment( (i-1)*nState+1 : (i-1)*nState+nState, nOdd*nState+nControl+1 ) = ...
            (dX*Ainv*D(:,i) - ...
            1/2*1/nSegment*dyn(1,C*B(:,i),uSegment) - ...
            Deltat/2*dxGradEven(:,1:nState,i)*dX*Ainv*B(:,i))*W(i);
    
    end
    
end

end