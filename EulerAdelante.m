% Euler hacia adelante

function [u,v,p,q] = EulerAdelante (h,u0,v0,p0,q0)
  
  pasos = 100/h;
  
	u = zeros(pasos+1,1) ;
	v = zeros(pasos+1,1) ;
	p = zeros(pasos+1,1) ;
	q = zeros(pasos+1,1) ;
	u(1) = u0 ;
	v(1) = v0 ;
	p(1) = p0 ;
	q(1) = q0 ;
  
  for i = 1 : pasos
    u(i+1) = u(i) + h*p(i) ;
    v(i+1) = v(i) + h*q(i) ;
    p(i+1) = p(i) + h*(-2)*v(i)*p(i)*q(i)/((u(i))^2+(v(i))^2+1) ;
    q(i+1) = q(i) + h*(-2)*u(i)*p(i)*q(i)/((u(i))^2+(v(i))^2+1) ;
  end  
  
end

