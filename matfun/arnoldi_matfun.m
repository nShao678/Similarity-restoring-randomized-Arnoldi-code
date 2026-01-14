function [y, hist, info] = arnoldi_matfun(A, q, f, m, para)
	% [y, hist, info] = arnoldi_matfun(A, q, m, para)
	% Computes f(A)q with Arnoldi, randomized Arnoldi, or oblique Arnoldi
	% Input:
	% A: matrix, n x n
	% q: vector, n x 1
	% f: function
	% m: iterations for Arnoldi
	% para.method: 'standard', 'rand', 'oblique'
	% para.tol: tolerance for solution
	% para.truesol: true solution for error comparison [optional]
	% para.sketch.d: sketching dimension
	% para.sketch.Omega: sketching matrix (function handle)
	% para.spacing: compute solution once every para.spacing iterations 
	% para.lastorth_method: 'chol' or 'lsqr' (method for solving least squares correction in oblique method)
	% para.lsqr_tol: tolerance for LSQR
	% Output:
	% y: solution, n x 1
	% hist: history of solutions, n x m [optional]
	% info: struct with extra information [optional]
	%	info.Q	basis
	%	info.H	Hessenberg matrix

	if nargin < 5
		para = [];
	end
	
	if ~isfield(para,'sketch')
		para.method = 'standard';
	else
		d = para.sketch.d;
		sketch = para.sketch.Omega;
	end
	if ~isfield(para, 'lastorth_method')
		para.lastorth_method = 'chol';
	end
	if ~isfield(para, 'spacing')
		para.spacing = 1;
	end
	if ~isfield(para, 'lsqr_tol')
		para.lsqr_tol = 1e-12;
	end
	if ~isfield(para, 'tol')
		para.tol = 0;		% never stop early
	end
	if ~isfield(para, 'truesol')
		get_trueerr = false;
	else
		get_trueerr = true;
	end

	n = length(q);
	Q = zeros(n, m+1);
	H = zeros(m+1, m);
	e1 = zeros(m, 1);	e1(1) = 1;
	if strcmp(para.method,'standard')
		nq_init = norm(q);
		q = q/nq_init;
		Q(:, 1) = q;
	else
		sQ = zeros(d, m);
		sq = sketch(q);
		nsq_init = norm(sq);
		q = q/nsq_init;
		Q(:, 1) = q;
		sQ(:,1) = sq/nsq_init;
	end

	if nargout >= 2
		nsols = ceil(m/para.spacing);
		hist.iter = [para.spacing*(1:nsols-1), m];
		hist.sol = zeros(n, nsols);
		hist.err = zeros(1, nsols);
		hist.err(1) = Inf;
	else
		hist.iter = m;
		hist.sol = zeros(n, 1);
		hist.err = Inf;
	end
	if nargout == 3
		info.ev = cell(nsols);
		info.evflag = zeros(1, nsols);
		info.evflagsmall = zeros(1, nsols);
		info.evflagcomplex = zeros(1, nsols);
	end
		
	solidx = 1;		% index in hist.iter of next solution that will be computed
	for j = 1:m
		if mod(j, 50) == 0
			fprintf("%d ", j);
		end
		% expand basis
		w = A(q);
        if strcmp(para.method,'standard')
			% cgs2
            h = Q(:,1:j)'*w;
            w = w-Q(:,1:j)*h;
            hh = Q(:,1:j)'*w;
            q = w-Q(:,1:j)*hh;
            h = h+hh;
            beta = norm(q);
        else
            sw = sketch(w);
			% rgs
            % h = sQ(:,1:j)\sw;
			h = sQ(:,1:j)'*sw;
			q = w - Q(:,1:j)*h;
			% rgs2 (reorth) -- unnecessary
			sq = sketch(q);
			hh = sQ(:,1:j)'*sq;
            q = q - Q(:,1:j)*hh;
			h = h + hh;

            sq = sketch(q);
            beta = norm(sq);
            sq = sq/beta;
            sQ(:, j+1) = sq;
        end
		q = q/beta;
        Q(:, j+1) = q;
        H(1:j, j) = h;
        H(j+1, j) = beta;


		if j == hist.iter(solidx)
			% compute approximate solution
			if strcmp(para.method,'standard')
				% orthogonal basis
				[V, D] = eig(H(1:j,1:j));	% H = V*D*V^{-1}
				fD = diag(f(diag(D)));
				hist.sol(:, solidx) = Q(:, 1:j) * (V * (fD * (V \ (nq_init*e1(1:j)))));
			elseif strcmp(para.method,'rand')
				% sketch-orthogonal basis
				[V, D] = eig(H(1:j,1:j));	% H = V*D*V^{-1}
				fD = diag(f(diag(D)));
				hist.sol(:, solidx) = Q(:, 1:j) * (V * (fD * (V \ (nsq_init*e1(1:j)))));
			else	% oblique Arnoldi
				% compute oblique Arnoldi decomposition
				[h, info_lsqr] = lastcol_orth(Q(:,1:j+1), para.lastorth_method, para.lsqr_tol);
				if (strcmp(para.lastorth_method, 'lsqr') && nargout == 3)
					info.lsqr_iter(solidx) = info_lsqr.iter;
					info.lsqr_res(solidx) = info_lsqr.lsvec(end);
				end
		
				Hnew = H(1:j,1:j);
				Hnew(1:j,j) = Hnew(1:j,j) + beta*h;
				if nargout == 3
					info.Hupdate{solidx} = beta*h;
				end
				% compute solution
				[V, D] = eig(Hnew);	% H = V*D*V^{-1}
				fD = diag(f(diag(D)));
				hist.sol(:, solidx) = Q(:, 1:j) * (V * (fD * (V \ (nsq_init*e1(1:j)))));
			end
			% Check convergence:
			if para.tol > 0
				if get_trueerr
					% true error
					hist.err(solidx) = norm(hist.sol(:, solidx) - para.truesol)./norm(para.truesol);
				else
					% consecutive difference
					if solidx == 1
						hist.err(solidx) = Inf;
					else
						hist.err(solidx) = norm(hist.sol(:,solidx) - hist.sol(:,solidx-1))./norm(hist.sol(:,solidx));
					end
				end
				% Compare with tolerance:
				if hist.err(solidx) < para.tol
					hist.err = hist.err(1:solidx);
					hist.sol = hist.sol(:, 1:solidx);
					hist.iter = hist.iter(:, 1:solidx);
					y = hist.sol(:, end);
					break;
				end
			end
			solidx = solidx + 1;
		end
		y = hist.sol(:, end);		% final solution
	end

	if nargout == 3
		info.H = H;
		info.Q = Q;
		% if strcmp(para.method, "oblique")
		% 	info.Hnew = Hnew;
		% end
		if ~strcmp(para.method, "standard")
			info.sQ = sQ;
		end
	end


end
