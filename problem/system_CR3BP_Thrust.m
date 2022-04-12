function [dx,dxGrad] = system_CR3BP_Thrust(t,x,u,param)

% - CR3BP + Thrust dynamics
[dx] = dynamics_CR3BP_Thrust(t,x,u,param);

if nargout == 2 
    % - analytic gradients
    dxGrad = jacobian_CR3BP_Thrust(t,x,u,param);   
end

end