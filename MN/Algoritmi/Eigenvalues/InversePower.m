% [USES] Gauss/G/G.m
% [USES] Algoritmi-ad-hoc/SST.m
% [USES] Algoritmi-ad-hoc/SIT.m
function [k, l , y] = InversePower(A, max_iter, tol, miu)
	% This function returns the smallest eigenvalue of the matrix and 
	% the characteristic vector
	n = length(A);
	% generate a random vector of length n
	y0 = rand(n, 1);
	% generate the I matrix of length n
	I = eye(n);
	% iterate from 1 to the maximum number of iterations
	for k = 1:max_iter
		% caculate z by solving the system using the inverse matrix method
		B = A - miu * I;
		z = G (B, y0);
		% calculate the eigenvector associated to the eigenvalue of the
		% current step
		y = z / norm(z);
		% calculate the eigenvalue at the current step
		l = y' * A * y;
		% to speed up the convergence, starting from the second step
		% miu's value becomes the value of the eigenvalue calculated
		% at the current step
		if k > 1
			miu = l;
		endif
		% calculate the relative variance(epsilon) at the current step
		eps = abs((y - y0) / y);
			% if epsilon is smaller then the tolerance we end the algorithm
		if eps < tol
			break;
		endif
		% y_k from the current step becomes y_(k - 1) for the next step
		y0 = y;
	endfor
end
