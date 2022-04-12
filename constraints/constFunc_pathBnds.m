function Fpathbnds = constFunc_pathBnds(z,transcribe,problem)

% - extract parameters
nState = problem.nState;
nSegment = transcribe.param.nSegment;
nOdd = transcribe.param.nOdd;

% - path end point state bounds
Fpathbnds = [ 
    z(1:nState);                                                    % initial state bound
    z(nState*(nOdd-1)*nSegment+1:nState*(nOdd-1)*nSegment+nState);  % final state bound 
    ];

end