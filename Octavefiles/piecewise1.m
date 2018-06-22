function [y] = piecewise1(x)
  y = 0;
  if x >= 1 && x <= 4
    y = -x + 8;
  elseif x > 4 && x <= 6
    y = 2*x - 4;
  endif
end function