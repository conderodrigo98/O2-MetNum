% lsode

function xdot = f (x, t)

  xdot = zeros (4,1);

  xdot(1) = x(3);
  xdot(2) = x(4);
  xdot(3) = (-2*x(2)*x(3)*x(4))/(x(1)^2 + x(2)^2 +1);
  xdot(4) = (-2*x(1)*x(3)*x(4))/(x(1)^2 + x(2)^2 +1);

endfunction
