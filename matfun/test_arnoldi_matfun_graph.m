
clear;

addpath('../utils/');

numtests = 1;
graphnames = ["big", "foldoc", "p2p-Gnutella30"];
times = zeros(length(graphnames), numtests, 5);
iters = zeros(length(graphnames), numtests, 5);
errors = zeros(length(graphnames), numtests, 5);
datastring = strings(length(graphnames), 1);
rng(0);

for jj = 1:length(graphnames)
	load(graphnames(jj));

	A = Problem.A;
	A = (A~=0);		% take sparsity pattern only
	A = diag(sum(A)) - A;
	n = size(A, 1);
	
	b = randn(n, 1);
	b = b/norm(b);
	
	A = @(x) A*x;

	m = 500;
	d = 1000;

	Omega = sparse_sign(d,n,8);
	para.sketch.Omega = @(x) Omega*x;
	para.sketch.d = d;

	f = @(z) sqrt(z);
	testlabel = 'sqrt';
	testidx = 1;

	% "true" solution
	% y_true = idct(f(A0).*dct(b));
	para_true.method = 'standard';
	para_true.spacing = 25;
	para_true.tol = 1e-8;
	[y_true, hist_true] = arnoldi_matfun(A, b, f, 1000, para_true);
	fprintf("\nFinal sol:\tError: %.2e\tIter: %d\n", hist_true.err(end), hist_true.iter(end));

	runs = 1;	% multiple runs to get average times
	for run = 1:runs
		fprintf('%d\n', run);
		para.spacing = 10;
		para.tol = 1e-6;
		para.truesol = y_true;
		para.method = 'standard';
		fprintf("%s:\t", para.method);
		tic;
		[y_std, hist_std] = arnoldi_matfun(A, b, f, m, para);
		times(jj, testidx, 1) = times(jj, testidx, 1) + toc/runs;
		iters(jj, testidx, 1) = hist_std.iter(end);
		errors(jj, testidx, 1) = hist_std.err(end);
		fprintf("\n");
		nerr_std{jj} = vecnorm(hist_std.sol - y_true)./norm(y_true);
		iter_std{jj} = hist_std.iter;

		para.method = 'rand';
		fprintf("%s:\t\t", para.method);
		tic;
		[y_rand, hist_rand] = arnoldi_matfun(A, b, f, m, para);
		times(jj, testidx, 2) = times(jj, testidx, 2) + toc/runs;
		iters(jj, testidx, 2) = hist_rand.iter(end);
		errors(jj, testidx, 2) = hist_rand.err(end);
		fprintf("\n");
		nerr_rand{jj} = vecnorm(hist_rand.sol - y_true)./norm(y_true);
		iter_rand{jj} = hist_rand.iter;

		para.method = 'oblique';
		para.lastorth_method = 'chol';
		fprintf("%s:\t", para.method);
		tic;
		[y_obl, hist_obl] = arnoldi_matfun(A, b, f, m, para);
		times(jj, testidx, 3) = times(jj, testidx, 3) + toc/runs;
		iters(jj, testidx, 3) = hist_obl.iter(end);
		errors(jj, testidx, 3) = hist_obl.err(end);
		fprintf("\n");
		nerr_obl{jj} = vecnorm(hist_obl.sol - y_true)./norm(y_true);
		iter_obl{jj} = hist_obl.iter;

		para.method = 'oblique';
		para.lastorth_method = 'lsqr';
		para.lsqr_tol = 1e-12;
		fprintf("%s:\t", para.method);
		tic;
		[y_obl2, hist_obl2] = arnoldi_matfun(A, b, f, m, para);
		times(jj, testidx, 4) = times(jj, testidx, 4) + toc/runs;
		iters(jj, testidx, 4) = hist_obl2.iter(end);
		errors(jj, testidx, 4) = hist_obl2.err(end);
		fprintf("\n");
		nerr_obl2{jj} = vecnorm(hist_obl2.sol - y_true)./norm(y_true);
		iter_obl2{jj} = hist_obl2.iter;

		para.method = 'oblique';
		para.lastorth_method = 'lsqr';
		para.lsqr_tol = 1e-1;
		fprintf("%s:\t", para.method);
		tic;
		[y_obl3, hist_obl3] = arnoldi_matfun(A, b, f, m, para);
		times(jj, testidx, 5) = times(jj, testidx, 5) + toc/runs;
		iters(jj, testidx, 5) = hist_obl3.iter(end);
		errors(jj, testidx, 5) = hist_obl3.err(end);
		fprintf("\n");
		nerr_obl3{jj} = vecnorm(hist_obl3.sol - y_true)./norm(y_true);
		iter_obl3{jj} = hist_obl3.iter;

	end
end

clear hist_true hist_std hist_rand hist_obl hist_obl2 hist_obl3;
clear y_true y_std y_rand y_obl y_obl2 y_obl3;
save("data_graph");

%% 
% Figures with error: 
for jj = 1:length(graphnames)
	leg = strings(0);
	figure;
	semilogy(iter_std{jj}, nerr_std{jj}, 'k-', 'LineWidth',2);	hold on;
	leg(end+1) = "Arnoldi";
	
	semilogy(iter_rand{jj}, nerr_rand{jj}, 'b-', 'LineWidth',2);	hold on;
	leg(end+1) = "randomized Arnoldi";
	
	semilogy(iter_obl{jj}, nerr_obl{jj}, 'r-', 'LineWidth',2);	hold on;
	leg(end+1) = "SRR-Arnoldi (chol)";
	
	semilogy(iter_obl2{jj}, nerr_obl2{jj}, 's--', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (lsqr $10^{-12}$)";
	
	semilogy(iter_obl3{jj}, nerr_obl3{jj}, 'mo-.', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (lsqr $10^{-1}$)";
	
	legend(leg, 'fontsize', 26, 'location', 'southwest', 'interpreter', 'latex', 'box', 'off');
	maxiter = max([max(iter_std{jj}), max(iter_rand{jj}), max(iter_obl{jj}), max(iter_obl2{jj}), max(iter_obl3{jj})]);
	xlim([-15, maxiter + 15]);
	ylim([1e-9, 1e-1]);
	title(graphnames(jj), 'fontsize',26);
	xlabel('iteration', 'fontsize',20);
	ylabel('relative error', 'fontsize',20);
	drawnow;

	ax = gca;
	ax.XAxis.FontSize = 20;
	ax.YAxis.FontSize = 20;

end


% Figure with times:
for ii = 1:numtests
	figure;
	T = squeeze(times(:, ii, :));
	b = bar(T);

	b(1).FaceColor = 'k';
	b(2).FaceColor = 'b';
	b(3).FaceColor = 'r';
	b(5).FaceColor = 'm';

	xticks(1:length(graphnames));
	xticklabels(graphnames);
	xlim([0.5, length(graphnames) + 0.5]);
	ylim([0, max(T(:))*1.5]);
	ylabel('time (s)', 'fontsize',20);
	title('Times', 'fontsize',20);
	legend(leg, 'FontSize', 20, 'Location', 'northwest', 'interpreter', 'latex', 'box','off');

	drawnow;

	ax = gca;
	ax.XAxis.FontSize = 20;
	ax.YAxis.FontSize = 20;

end

