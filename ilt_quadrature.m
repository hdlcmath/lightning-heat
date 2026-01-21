function [qW,qP] = ilt_quadrature(t,M)
% Given a value of t, output a vector of weights and M points for inverse laplace transform

th = (.5:M-.5)*pi/M; % Vector of theta to go in parametrization
sg = -0.6122; mu = 0.5017; bt = 0.6407; nu = 0.2645;
z=@(th,t) 2*N./t.*(sg + mu*th.*cot(bt*th)+nu*1i*th);
dz=@(th,t) 2*N./t.*(mu*cot(bt*th)...
		-mu*bt*th.*csc(bt*th).^2+nu*1i);
	
qP = z(th,t);
qW = exp(z(th,t)*t).*dz(th,t);

end
