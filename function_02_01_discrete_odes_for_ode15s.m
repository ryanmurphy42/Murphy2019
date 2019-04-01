function ut=function_02_01_discrete_odes_for_ode15s(t,u)
%function_02_01_discrete_odes_for_ode15s - ode for each spring boundary

% Problem parameters
  global ncall
  global xu
  global xl
  global eta_spring
  global k_spring
  global a_spring

 
  n=length(u);
  ut = zeros(1,n);

  for i=1:n
    if (i==1)    
		ut(i) = 0; %The first position is fixed.
    elseif (i==n)
         ut(i) = 0; %The final position is fixed.
    else
		ut(i) = -(k_spring(i-1)/eta_spring)*(u(i)-u(i-1) -a_spring(i-1)) + (k_spring(i)/eta_spring)*(u(i+1)-u(i) -a_spring(i));
    end
  end
  ut = ut';
      
     
% Increment calls to pde_1
  ncall=ncall+1;
  
end