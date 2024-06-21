function [ t, y ] = odeEuler( ode, tspan, y0, options, varargin )
%ODEEULER
%

% This is an extremely simple way of solving an ODE
% Uses formula $y_{n+1} = y_n + \delta t * f(t_n, y_n)$

% Evan's Test Code
% function dydt = simpleODE(t, y)
%     dydt = -y + t;
% end
% options = struct('TimeStep', 0.1, 'UseOutputFcn', false);
% [t, y] = odeEuler(@simpleODE, [0,10], 1, options)
% EVAN: I added these because I was testing out this function
%   and it kept giving me errors.
if ~isfield(options, 'OutputFcn')
    options.OutputFcn = [];
end
if ~isfield(options, 'OutputFcnEvalInterval')
    options.OutputFcnEvalInterval = inf; % No interval if not specified
end
if ~isfield(options, 'Precompute')
    options.Precompute = [];
end

dt = options.TimeStep;
UseOutputFcn = options.UseOutputFcn;
OutputFcn = options.OutputFcn;
OutputFcnEvalInterval = options.OutputFcnEvalInterval;

t0 = tspan(1);
tf = tspan(2);

[numSteps, numSamples, t] = odeEuler_tspan(tspan, dt);

n = size(y0, 1);

y = zeros(numSamples,n);

y(1, :) = y0;

if isempty(options.Precompute)
    precomputed_args = [];
else
    precomputed_args = options.Precompute(t);
    varargin{end+1} = [];
end

OutputFcnEvalTime = t0;
if UseOutputFcn
    OutputFcn(t0,y0,'init');
end

for i = 1:numSteps
    
    if ~isempty(precomputed_args)
        varargin{end} = precomputed_args(i,:);
    end

    try
        y(i+1,:) = y(i,:) + dt * ode(t(i+1), y(i,:)', varargin{:})';
    catch exception
        disp("t is a scalar"); disp(isscalar(t(i+1))); disp(t(i+1));
        throw(exception)
    end

    if UseOutputFcn
        if t(i+1) - OutputFcnEvalTime >= OutputFcnEvalInterval
            status = OutputFcn(t(i+1),y(i+1,:)',[]);
            if status == true
                break
            end
        end
    end
    
end

if UseOutputFcn
    OutputFcn([],[],'done');
end

end