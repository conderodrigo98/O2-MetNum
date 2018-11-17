% Heun

function [u,v,p,q] = Heun (h, u0, v0, p0, q0)
	%condiciones iniciales
	pasos = 100/h;
	u = zeros(pasos+1,1);
	v = zeros(pasos+1,1);
	p = zeros(pasos+1,1);
	q = zeros(pasos+1,1);
	u(1) = u0;
	v(1) = v0;
	p(1) = p0;
	q(1) = q0;
	%iteración
  for i = 2 : pasos+1
		%predicción
		prediccionU = u(i-1) + h*p(i-1);
		prediccionV = v(i-1) + h*q(i-1);
		prediccionP = p(i-1) + h*(-2*v(i-1)*p(i-1)*q(i-1)/(u(i-1)^2*v(i-1)^2+1)); 
		prediccionQ = q(i-1) + h*(-2*u(i-1)*p(i-1)*q(i-1)/(u(i-1)^2*v(i-1)^2+1));
		%corrección
		u(i) = u(i-1) + (h/2)*(p(i-1) + prediccionP);
		v(i) = v(i-1) + (h/2)*(q(i-1) + prediccionQ);
		p(i) = p(i-1) + (h/2)*((-2*v(i-1)*p(i-1)*q(i-1)/(u(i-1)^2*v(i-1)^2+1) + (-2*prediccionV*prediccionP*prediccionQ/(prediccionU^2*prediccionV^2+1))));
		q(i) = q(i-1) + (h/2)*((-2*u(i-1)*p(i-1)*q(i-1)/(u(i-1)^2*v(i-1)^2+1) + (-2*prediccionU*prediccionP*prediccionQ/(prediccionU^2*prediccionV^2+1))));
	end
endfunction

