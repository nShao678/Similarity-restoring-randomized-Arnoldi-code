function [h, info] = lastcol_orth(Q, method, lsqr_tol)
	% [h, info] = lastcol_orth(Q, method, lsqr_tol)
	% method: "lsqr" or "chol"

	if nargin < 3
		lsqr_tol = 1e-12;
	end

	m = size(Q, 2);
	if strcmp(method, 'chol')
		M = Q(:, 1:m-1)'*Q;
		R = chol(M(:,1:m-1));
		h = R\((R')\M(:,m));
		if nargout == 2
			info = [];
		end
	elseif strcmp(method, 'lsqr')
		maxit = 100;
		[h, ~, ~, iter, ~, lsvec] = lsqr(Q(:, 1:m-1), Q(:, m), lsqr_tol, maxit);
		if nargout == 2
			info.iter = iter;
			info.lsvec = lsvec;
		end
	end

end
