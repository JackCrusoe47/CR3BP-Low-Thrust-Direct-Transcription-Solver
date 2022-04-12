function g = gradFunc(z,auxdata)
% -- optimzation gradient

problem = auxdata.problem;
transcribe = auxdata.transcribe;

switch problem.objType
    case 'propOptimal'
        g = gradFunc_propOptimal(z,transcribe,problem);
    case 'timeOptimal'
        g = gradFunc_timeOptimal(z,transcribe,problem);
    case 'thrustOptimal'
        g = gradFunc_thrustOptimal(z,transcribe,problem);
    case 'massTarget'
        g = gradFunc_massTarget(z,transcribe,problem);
    case 'stateTarget'
        g = gradFunc_stateTarget(z,transcribe,problem);
    case 'feasible'
        g = gradFunc_feasible(z,transcribe,problem);
    otherwise
        g = gradFunc_feasible(z,transcribe,problem);
end


end