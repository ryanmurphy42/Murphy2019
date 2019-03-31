function [x,isterm,dir] = function_04_eventfun(t,u)

 global ncall
  global xu
  global xl
  global eta_spring
  global k_spring
  global a_spring

 
  n=length(u);
  ut = zeros(1,n);
%event condition to stop when reaches steady state
for i=1:n
    if (i==1)
        ut(i) = 0; %The first position is fixed.
    elseif (i==n)
        ut(i) = 0; %The final position is fixed.
    else
        ut(i) = -(k_spring(i-1)/eta_spring)*(u(i)-u(i-1) -a_spring(i-1)) + (k_spring(i)/eta_spring)*(u(i+1)-u(i) -a_spring(i)); %Linear force law for internal spring boundaries
    end
end

dy = sum(abs(ut));
x = norm(dy) - 1e-10;
isterm = 1;
dir = 0;
end