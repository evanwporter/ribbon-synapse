function [Jconstant,Jfcn,Jargs,Joptions] = ...
    odejacobian(fcnHandlesUsed,ode,t,y0,options,extras)
%ODEJACOBIAN  Helper function for the Jacobian function in ODE solvers
%    ODEJACOBIAN determines whether the Jacobian is constant and if so,
%    returns its value as Jfcn. If an analytical Jacobian is available from
%    a function, ODEJACOBIAN initializes Jfcn and creates a cell array of
%    additional input arguments. For numerical Jacobian, ODEJACOBIAN tries to
%    extract JPattern and sets JOPTIONS for use with ODENUMJAC.
%
%   See also ODE15S, ODE23S, ODE23T, ODE23TB, ODENUMJAC.

%   Jacek Kierzenka
%   Copyright 1984-2020 The MathWorks, Inc.

Jconstant = strcmp(odeget(options,'JConstant','off'),'on');
Jfcn = [];
Jargs = {};
Joptions = [];   

Janalytic = false;

if fcnHandlesUsed
  Jfcn = odeget(options,'Jacobian',[]);
  if ~isempty(Jfcn)
    if isnumeric(Jfcn)
      Jconstant = true;
    else
      Janalytic = true;
      Jargs = extras;
    end
  end  
else  % ode-file used  
  joption = odeget(options,'Jacobian','off');  
  switch lower(joption)
    case 'on'    % ode(t,y,'jacobian',p1,p2...)
      Janalytic = true;
      Jfcn = ode;
      Jargs = [{'jacobian'} extras];  
    case 'off'   % use odenumjac
    otherwise
      error(message('MATLAB:odejacobian:InvalidJOption', joption));
  end    
end  

if ~Janalytic   % odenumjac will be used
  Joptions.diffvar  = 2;       % df(t,y)/dy
  Joptions.vectvars = [];  
  vectorized = strcmp(odeget(options,'Vectorized','off'),'on');  
  if vectorized
    Joptions.vectvars = 2;     % f(t,[y1,y2]) = [f(t,y1), f(t,y2)] 
  end
  
  atol = odeget(options,'AbsTol',1e-6);
  % Joptions.thresh = zeros(size(y0))+ atol(:);
  Joptions.thresh = zeros(size(y0))+ atol(t);
  Joptions.fac  = [];
  
  if fcnHandlesUsed  
    jpattern = odeget(options,'JPattern',[]); 
  else  % ode-file used
    jp_option = odeget(options,'JPattern','off'); 
    switch lower(jp_option)
      case 'on'
        jpattern = feval(ode,[],[],'jpattern',extras{:}); 
      case 'off'  % no pattern provided
        jpattern = [];
      otherwise
        error(message('MATLAB:odejacobian:InvalidJpOption', jp_option));        
    end          
  end  
  if ~isempty(jpattern)
    Joptions.pattern = jpattern;    
    Joptions.g = colgroup(jpattern);
  end  
end
    
