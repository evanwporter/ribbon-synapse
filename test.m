% function dydt = simpleODE(t, y)
%     dydt = -y + t;
% end

simpleODE = @(t, y) -y + t;

options = struct('TimeStep', 0.1, 'UseOutputFcn', false);

[t, y] = odeEuler(@simpleODE, [0,10], 1, options)
