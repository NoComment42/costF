function [y] = piecewise2(x)
  y = 0;
  if x >= 1 && x <= 4
    y = (1+1/3)*x + 2/3;
  elseif x > 3 && x <= 6
    y = -2.5*x + 16;
  endif
end function