function [u,matvec,hist,info,flag] = KS(A,q,numeval,para)

%para.dMax: The maximum allowed dimension of Krylov space
%para.dMin: The number of Ritz vectors kept after restarting
%para.type: 'largestabs' etc, same as eigs !! If using 'smallestabs', A(x)
%should be A\x
%para.getinfo: if true, also output info (more expensive, not good for testing runtime!)
%numeval: compute the numeval eigenvalues with the largest magnitude
%flag = 0: KS does not converge
if nargin==3
    para = [];
end

if ~isfield(para,'dMin')
    para.dMin = 2*numeval; % restarting cycle
end
if ~isfield(para,'dMax')
    para.dMax = 4*numeval; % restarting cycle
end
if ~isfield(para,'matvecMax')
    para.matvecMax = 2000; % restarting cycle
end
if ~isfield(para,'tolr')
    para.tolr = 1e-6;   % SC: norm(R,fro)/nA<tolr
end
if ~isfield(para,'type')
    para.type = 'largestabs';   % same as 'eigs'
end
if ~isfield(para,'sketch')
    para.method = 'KS';
else
    d = para.sketch.d;
    sketch = para.sketch.Omega;
end
if ~isfield(para, 'getinfo')
	para.getinfo = false;
end
if ~isfield(para, 'lastorth_method')
	para.lastorth_method = 'chol';
end
if ~isfield(para, 'lsqr_tol')
	para.lsqr_tol = 1e-12;
end
if ~isfield(para, 'ks_reorth')
	para.ks_reorth = true;		% if true, use cgs2 in KS; else, use cgs
end


[n,~] = size(q);
Q = zeros(n, para.dMax);
H = zeros(para.dMax, para.dMax-1);
if strcmp(para.method,'KS')
    q = q/norm(q);
    Q(:, 1) = q;
else
    sQ = zeros(d,para.dMax);
    sq = sketch(q);
    nsq = norm(sq);
    q = q/nsq;
    Q(:, 1) = q;
    sQ(:,1) = sq/nsq;
end
j0 = 1;
hist.eval = zeros(numeval,ceil(para.matvecMax/(para.dMax-para.dMin)));
hist.res = zeros(numeval,ceil(para.matvecMax/(para.dMax-para.dMin)));
matvec = 0;
flag = [];
info = [];

for iter = 1:para.matvecMax
    for j = j0:para.dMax-1
        w = A(q);
        if strcmp(para.method,'KS')
			% cgs
            h = Q(:,1:j)'*w;
            q = w-Q(:,1:j)*h;
			if para.ks_reorth
				% cgs2
				hh = Q(:,1:j)'*q;
				q = q-Q(:,1:j)*hh;
				h = h+hh;
			end
            beta = norm(q);
        else
            sw = sketch(w);
            h = sQ(:,1:j)\sw;
            q = w-Q(:,1:j)*h;
            sq = sketch(q);
            beta = norm(sq);
            sq = sq/beta;
            sQ(:, j+1) = sq;
        end
        q = q/beta;
        Q(:, j+1) = q;
        H(1:j,j) = h;
        H(j+1,j) = beta;
    end

	if para.getinfo
		info.condQ(iter) = cond(Q(:, 1:para.dMax));
    end

    if strcmp(para.method,'OKS')
		% Least squares solve: 
		if para.getinfo
			[h, orthinfo] = lastcol_orth(Q, para.lastorth_method, para.lsqr_tol);
			if strcmp(para.lastorth_method, 'lsqr')
				info.lsqr_iter(iter) = orthinfo.iter;
				info.lsqr_res(iter) = orthinfo.lsvec(end);
			end
		else
			h = lastcol_orth(Q, para.lastorth_method, para.lsqr_tol);
		end
        % M = Q(:,1:para.dMax-1)'*Q;
        % R = chol(M(:,1:para.dMax-1));
        % h = R\((R')\M(:,para.dMax));
%         h = M(:,1:para.dMax-1)\M(:,para.dMax);

		% Update Q and H:
        q = q-Q(:,1:para.dMax-1)*h;
        sq = sq-sQ(:,1:para.dMax-1)*h;
        beta0 = norm(sq);
        q = q/beta0;
        sq = sq/beta0;
        Q(:, para.dMax) = q;
        sQ(:, para.dMax) = sq;
        H(1:para.dMax-1,para.dMax-1) = H(1:para.dMax-1,para.dMax-1)+beta*h;
        beta = beta*beta0;
        H(para.dMax,para.dMax-1) = beta;

		if para.getinfo
			info.condQ_after(iter) = cond(Q(:, 1:para.dMax));
		end
    end

    matvec = matvec+para.dMax-j0;

    [V,S] = schur(H(1:para.dMax-1,1:para.dMax-1));
    [VV,eval] = eig(S,'vector');
    VV = V*VV;
    res = abs(beta*VV(para.dMax-1,:))*norm(q);
    idx = whichEigenvalues(eval, para.type);
    eval = eval(idx);
    res = res(idx);
    hist.eval(:,iter) = eval(1:numeval);
    hist.res(:,iter) = res(1:numeval);
    
%     isNotConverged = ~(res(1:numeval)' < para.tolr*max(eps^(2/3), abs(eval(1:numeval))));
    isNotConverged = ~(res(1:numeval)' < para.tolr);
    nconv = nnz(~isNotConverged);

    if nconv >= numeval || matvec>=para.matvecMax
        u = Q(:,1:para.dMax-1)*VV(:,idx(1:numeval));
        if matvec>=para.matvecMax
            warning(['May not converge in ',num2str(iter),'cycles'])
            flag = 0;
        else
            flag = 1;
        end
        if para.hist == 1
            hist.eval = hist.eval(:,1:iter);
            hist.res = hist.res(:,1:iter);
        end
        break
    end

    eval = ordeig(S);
    idx = whichEigenvalues(eval,para.type);
    idx = idx(1:para.dMin);
    select = false(1,para.dMax-1);
    select(idx) = true;

    dMin = para.dMin;
    if isreal(H)
        for i = idx'
            if i < para.dMax-1 && S(i+1,i) ~= 0 && ~select(i+1)
                select(i+1) = true;
                dMin = dMin+1;
            end
            if i > 1 && S(i, i-1) ~= 0 && ~select(i-1)
                select(i-1) = true;
                dMin = dMin+1;
            end
        end
    end

    [V,S] = ordschur(V,S,select);
    if ~strcmp(para.method,'KS')
        sQ(:,1:dMin) = sQ(:,1:para.dMax-1)*V(:,1:dMin);
        sQ(:,dMin+1) = sq;
    end
    Q(:,1:dMin) = Q(:,1:para.dMax-1)*V(:,1:dMin);
    Q(:,dMin+1) = q;

    H = zeros(para.dMax, para.dMax-1);
    H(1:dMin+1,1:dMin) = [S(1:dMin,1:dMin);beta*V(para.dMax-1,1:dMin)];
    j0 = dMin+1;
end


end

function ind = whichEigenvalues(d, method)

switch method
    case 'largestabs'
        [~, ind] = sort(abs(d), 'descend');
    case 'largestreal'
        [~, ind] = sort(real(d), 'descend');
    case 'smallestreal'
        [~, ind] = sort(real(d), 'ascend');
    case 'largestimag'
        [~, ind] = sort(imag(d), 'descend');
    case 'smallestimag'
        [~, ind] = sort(imag(d), 'ascend');
    case 'bothendsreal'
        [~, ind] = sort(real(d), 'descend');
        ind2 = [ind, flip(ind)]';
        ind2 = ind2(:);
        ind = ind2(1:size(d,1));
    case 'bothendsimag'
        [~, ind] = sort(imag(d), 'descend');
        ind2 = [ind, flip(ind)]';
        ind2 = ind2(:);
        ind = ind2(1:size(d,1));
    case 'smallestabs'
        [~,ind] = sort(abs(d), 'ascend');
    case 'smallestimagabs'
        [~,ind] = sort(abs(imag(d)), 'ascend');
end

end
