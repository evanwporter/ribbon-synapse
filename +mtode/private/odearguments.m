function [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, odeFcn, ...
          options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, ...
          dataType ] =   ...
    odearguments(FcnHandlesUsed, solver, ode, tspan, y0, options, extras)
%ODEARGUMENTS  Helper function that processes arguments for all ODE solvers.
%
%   See also ODE113, ODE15I, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Mike Karr, Jacek Kierzenka
%   Copyright 1984-2020 The MathWorks, Inc.

if strcmp(solver,'ode15i')
  FcnHandlesUsed = true;   % no MATLAB v. 5 legacy for ODE15I
end  

if FcnHandlesUsed  % function handles used
  if isempty(tspan) || isempty(y0) 
    error(message('MATLAB:odearguments:TspanOrY0NotSupplied', solver));
  end      
  if length(tspan) < 2
    error(message('MATLAB:odearguments:SizeTspan', solver));
  end  
  htspan = abs(tspan(2) - tspan(1));  
  tspan = tspan(:);
  ntspan = length(tspan);
  t0 = tspan(1);  
  next = 2;       % next entry in tspan
  tfinal = tspan(end);     
  args = extras;                 % use f(t,y,p1,p2...) 

else  % ode-file used   (ignored when solver == ODE15I)
  % Get default tspan and y0 from the function if none are specified.
  if isempty(tspan) || isempty(y0) 
    if exist(ode)==2 && ( nargout(ode)<3 && nargout(ode)~=-1 )  %#ok<EXIST>
      error(message('MATLAB:odearguments:NoDefaultParams', funstring( ode ), solver, funstring( ode )));      
    end
    [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
    if isempty(tspan)
      tspan = def_tspan;
    end
    if isempty(y0)
      y0 = def_y0;
    end
    options = odeset(def_options,options);
  end  
  tspan = tspan(:);
  ntspan = length(tspan);
  if ntspan == 1    % Integrate from 0 to tspan   
    t0 = 0;          
    next = 1;       % Next entry in tspan.
  else              
    t0 = tspan(1);  
    next = 2;       % next entry in tspan
  end
  htspan = abs(tspan(next) - t0);
  tfinal = tspan(end);   
  
  % The input arguments of f determine the args to use to evaluate f.
  if (exist(ode)==2) %#ok<EXIST>
    if (nargin(ode) == 2)           
      args = {};                   % f(t,y)
    else
      args = [{''} extras];        % f(t,y,'',p1,p2...)
    end
  else  % MEX-files, etc.
    try 
      args = [{''} extras];        % try f(t,y,'',p1,p2...)     
      feval(ode,tspan(1),y0(:),args{:});   
    catch
      args = {};                   % use f(t,y) only
    end
  end
end

y0 = y0(:);
neq = length(y0);

% Test that tspan is internally consistent.
if any(isnan(tspan))
  error(message('MATLAB:odearguments:TspanNaNValues'));
end
if t0 == tfinal
  error(message('MATLAB:odearguments:TspanEndpointsNotDistinct'));
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
  error(message('MATLAB:odearguments:TspanNotMonotonic'));
end

f0 = feval(ode,t0,y0,args{:});   % ODE15I sets args{1} to yp0.
[m,n] = size(f0);
if n > 1
  error(message('MATLAB:odearguments:FoMustReturnCol', funstring( ode )));
elseif m ~= neq
    error(message('MATLAB:odearguments:SizeIC', funstring( ode ), m, neq, funstring( ode )));
end

% Determine the dominant data type
classT0 = class(t0);
classY0 = class(y0);
classF0 = class(f0);
if strcmp(solver,'ode15i')  
  classYP0 = class(args{1});  % ODE15I sets args{1} to yp0.
  dataType = superiorfloat(t0,y0,args{1},f0);

  if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
        strcmp(classF0,dataType) && strcmp(classYP0,dataType))
    input1 = '''t0'', ''y0'', ''yp0''';
    input2 = '''f(t0,y0,yp0)''';
    warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
  end    
else  
  dataType = superiorfloat(t0,y0,f0);
  
  if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
        strcmp(classF0,dataType))
    input1 = '''t0'', ''y0'''; 
    input2 = '''f(t0,y0)''';
    warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
  end        
end

% Get the error control options, and set defaults.
rtol = odeget(options,'RelTol',1e-3);
if (length(rtol) ~= 1) || (rtol <= 0)
  error(message('MATLAB:odearguments:RelTolNotPosScalar'));
end
if rtol < 100 * eps(dataType) 
  rtol = 100 * eps(dataType);
  warning(message('MATLAB:odearguments:RelTolIncrease', sprintf( '%g', rtol )))
end
atol = odeget(options,'AbsTol',1e-6);
if isa(atol, 'function_handle')
  normcontrol = strcmp(odeget(options,'NormControl','off'),'on');
  if normcontrol
    atol_try = atol(0);
    if length(atol_try) ~= 1
      error(message('MATLAB:odearguments:NonScalarAbsTol'));
    end
    normy = norm(y0);
  else
    atol_try = atol(0);
    if (length(atol_try) ~= 1) && (length(atol_try) ~= neq)
      error(message('MATLAB:odearguments:SizeAbsTol', funstring( ode ), neq)); 
    end
    assert(size(atol_try,1) == numel(atol_try), 'atol function must return a column vector.')
    normy = [];
  end
  threshold = @(t) atol(t) / rtol;
else
  if any(atol <= 0)
    error(message('MATLAB:odearguments:AbsTolNotPos'));
  end
  normcontrol = strcmp(odeget(options,'NormControl','off'),'on');
  if normcontrol
    if length(atol) ~= 1
      error(message('MATLAB:odearguments:NonScalarAbsTol'));
    end
    normy = norm(y0);
  else
    if (length(atol) ~= 1) && (length(atol) ~= neq)
      error(message('MATLAB:odearguments:SizeAbsTol', funstring( ode ), neq)); 
    end
    atol = atol(:);
    normy = [];
  end
  threshold = atol / rtol;
end
% By default, hmax is 1/10 of the interval.
safehmax = 16.0*eps(dataType)*max(abs(t0),abs(tfinal));  % 'inf' for tfinal = inf
defaulthmax = max(0.1*(abs(tfinal-t0)), safehmax);
hmax = min(abs(tfinal-t0), abs(odeget(options,'MaxStep',defaulthmax)));
if hmax <= 0
  error(message('MATLAB:odearguments:MaxStepLEzero'));
end
htry = abs(odeget(options,'InitialStep',[]));
if ~isempty(htry) && (htry <= 0)
  error(message('MATLAB:odearguments:InitialStepLEzero'));
end

odeFcn = ode;
