function f = objFunc(z,auxdata)
% -- optimzation objective

problem = auxdata.problem;
transcribe = auxdata.transcribe;

switch problem.objType
    case 'propOptimal'
        f = objFunc_propOptimal(z,transcribe,problem);
    case 'timeOptimal'
        f = objFunc_timeOptimal(z,transcribe,problem);
    case 'thrustOptimal'
        f = objFunc_thurstOptimal(z,transcribe,problem);
    case 'massTarget'
        f = objFunc_massTarget(z,transcribe,problem);
    case 'stateTarget'
        f = objFunc_stateTarget(z,transcribe,problem);
    case 'feasible'
        f = objFunc_feasible(z,transcribe,problem);
    otherwise
        f = objFunc_feasible(z,transcribe,problem);
end


end