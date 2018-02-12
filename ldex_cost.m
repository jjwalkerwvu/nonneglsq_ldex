% ldex_cost.m

% this is a function that computes the cost and the gradient, in order to help 
% find the ldex response function.

% the way I have ena_arr and ldex_current saved, they need to be transposed.
function [J,grad]=ldex_cost(theta_coeffs,ena_arr,ldex_current,lambda);

  ldex_current=ldex_current';
  ena_arr=ena_arr';
  
 
  m=length(ldex_current);
  
  % remember: there is no theta_0 term here! (No offset)
  J=1/2/m*(ena_arr*theta_coeffs-ldex_current)'*(ena_arr*theta_coeffs-ldex_current)+...
    lambda/2/m*(theta_coeffs(1:end)'*theta_coeffs(1:end));
  
  % compute the gradient:
  
  % remember: there is no theta_0 term here! (No offset)
  grad=1/m*(ena_arr*theta_coeffs-ldex_current)'*ena_arr+(lambda/m*theta_coeffs(1:end))';
  % need to perform transpose, because the result of the above is a row vector 
  grad=grad';