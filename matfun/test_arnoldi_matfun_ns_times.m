
clear;

addpath('../utils/');

nn = [10 20 40 80]*1e3;

numtests = 3;
times = zeros(length(nn), numtests, 5);
iters = zeros(length(nn), numtests, 5);
errors = zeros(length(nn), numtests, 5);

coeffs = [0.37 0.32 0.30 0.29];

rng(1)

for jj = 1:length(nn)
	n = nn(jj);
	fprintf('\n%d\n', n);

	nclusters = 4;
	centers = logspace(0, 3, nclusters);
	variances = centers/100;
	clustersize = floor(n/nclusters);
	A0 = [];
	for j = 1:nclusters
		mu = centers(j);
		sigma = variances(j);
		temp = mu + randn(clustersize, 1)*sqrt(sigma);
		A0 = [A0; temp];
		if j == nclusters
			temp2 = mu + randn(n - nclusters*clustersize, 1)*sqrt(sigma);
			A0 = [A0; temp2];
		end
	end
	coeff = coeffs(jj);
	A1 = coeff*randn(n-1,1);	% non-symmetric

	A = @(x) fun(x, A0, A1);

	b = randn(n, 1);
	m = 1000;
	d = 1500;

	Omega = sparse_sign(d,n,8);
	para.sketch.Omega = @(x) Omega*x;
	para.sketch.d = d;

	for testidx = 1:numtests
		fprintf('Test: %d\n', testidx);
		switch testidx
			case 1
				f = @(z) sqrt(z);
				testlabel = 'sqrt';
				testlabels{1} = 'sqrt';
			case 2
				f = @(z) 1./sqrt(z);
				testlabel = 'invsqrt';
				testlabels{2} = 'invsqrt';
			case 3
				f = @(z) log(z);
				testlabel = 'log';
				testlabels{3} = 'log';
		end

		% "true" solution
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
			nerr_std = vecnorm(hist_std.sol - y_true)./norm(y_true);

			para.method = 'rand';
			fprintf("%s:\t\t", para.method);
			tic;
			[y_rand, hist_rand] = arnoldi_matfun(A, b, f, m, para);
			times(jj, testidx, 2) = times(jj, testidx, 2) + toc/runs;
			iters(jj, testidx, 2) = hist_rand.iter(end);
			errors(jj, testidx, 2) = hist_rand.err(end);
			fprintf("\n");
			nerr_rand = vecnorm(hist_rand.sol - y_true)./norm(y_true);

			para.method = 'oblique';
			para.lastorth_method = 'chol';
			fprintf("%s:\t", para.method);
			tic;
			[y_obl, hist_obl] = arnoldi_matfun(A, b, f, m, para);
			times(jj, testidx, 3) = times(jj, testidx, 3) + toc/runs;
			iters(jj, testidx, 3) = hist_obl.iter(end);
			errors(jj, testidx, 3) = hist_obl.err(end);
			fprintf("\n");
			nerr_obl = vecnorm(hist_obl.sol - y_true)./norm(y_true);

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
			nerr_obl2 = vecnorm(hist_obl2.sol - y_true)./norm(y_true);

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
			nerr_obl3 = vecnorm(hist_obl3.sol - y_true)./norm(y_true);
		end
	end
end

clear hist_true hist_std hist_rand hist_obl hist_obl2 hist_obl3;
clear y_true y_std y_rand y_obl y_obl2 y_obl3;
save("data_ns_times");		% save workspace

%%
leg = strings(0);
% figure;
% semilogy(hist_std.iter, nerr_std, 'k-', 'LineWidth',2);	hold on;
leg(end+1) = "Arnoldi";

% semilogy(hist_rand.iter, nerr_rand, 'b-', 'LineWidth',2);	hold on;
leg(end+1) = "randomized Arnoldi";

% semilogy(hist_obl.iter, nerr_obl, 'r-', 'LineWidth',2);	hold on;
leg(end+1) = "SRR-Arnoldi (chol)";

% semilogy(hist_obl2.iter, nerr_obl2, 's--', 'LineWidth',1.5);	hold on;
leg(end+1) = "SRR-Arnoldi (lsqr $10^{-12}$)";

% semilogy(hist_obl3.iter, nerr_obl3, 'mo-.', 'LineWidth',1.5);	hold on;
leg(end+1) = "SRR-Arnoldi (lsqr $10^{-1}$)";


% Figure with times:
for ii = 1:numtests
	% Barplot
	T = squeeze(times(:, ii, :));
	figure;
	b = bar(T);
	b(1).FaceColor = 'k';
	b(2).FaceColor = 'b';
	b(3).FaceColor = 'r';
	% b(4).FaceColor = [0.5 0.5 0.5];
	b(5).FaceColor = 'm';

	xticks(1:length(nn));
	xticklabels(nn);
	xlim([0.5, length(nn) + 0.5]);
	ylim([0, max(T(:))*1.45]);
	ylabel('time (s)','fontsize',20);
	xlabel('$n$', 'interpreter', 'latex','fontsize',20);
	titlestring = sprintf("Times - %s", testlabels{ii});
	title(titlestring,'FontSize',26);

	legend(leg, 'FontSize', 22, 'Location', 'northwest', 'interpreter', 'latex','box','off');

	drawnow;

	ax = gca;
	ax.XAxis.FontSize = 20;
	ax.YAxis.FontSize = 20;

end

% Figure with iters:
for ii = 1:numtests
	% Barplot
	I = squeeze(iters(:, ii, :));
	figure;
	b = bar(I);
	b(1).FaceColor = 'k';
	b(2).FaceColor = 'b';
	b(3).FaceColor = 'r';
	% b(4).FaceColor = [0.5 0.5 0.5];
	b(5).FaceColor = 'm';

	xticks(1:length(nn));
	xticklabels(nn);
	xlim([0.5, length(nn) + 0.5]);
	ylim([0, max(I(:))*1.7]);
	xlabel('$n$', 'interpreter', 'latex','fontsize',20);
	ylabel('iterations','fontsize',20);
	titlestring = sprintf("Iterations - %s", testlabels{ii});
	title(titlestring,'FontSize',26);
	legend(leg, 'FontSize', 22, 'Location', 'northwest', 'interpreter', 'latex','box','off');

	drawnow;

	ax = gca;
	ax.XAxis.FontSize = 20;
	ax.YAxis.FontSize = 20;

end

function y = fun(x,A0,A1)
	y = dct(x);
	y = A0.*y+[A1.*y(1:end-1);0];
	y = idct(y);
	end
