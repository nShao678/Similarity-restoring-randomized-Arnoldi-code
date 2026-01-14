clear;

addpath('../utils/');

n = 1e4;

rng(1)

nclusters = 4;
centers = logspace(0, 3, nclusters);
variances = centers/10;
clustersize = n/nclusters;
A0 = [];
for j = 1:nclusters
	mu = centers(j);
	sigma = variances(j);
	temp = mu + randn(clustersize, 1)*sqrt(sigma);
	A0 = [A0; temp];
end
A1 = zeros(n-1,1);	% symmetric

A = @(x) fun(x, A0, A1);

b = randn(n, 1);
m = 300;
d = 2*m;


Omega = sparse_sign(d,n,8);
para.sketch.Omega = @(x) Omega*x;
para.sketch.d = d;

numtests = 3;
times = zeros(numtests, 5);
testlabels = ["sqrt", "invsqrt", "log"];


for testidx = 1:numtests
	disp(testidx)
	switch testidx
		case 1
			f = @(z) sqrt(z);
			testlabel = 'sqrt';
		case 2
			f = @(z) 1./sqrt(z);
			testlabel = 'invsqrt';
		case 3
			f = @(z) log(z);
			testlabel = 'log';
	end

	% true solution
	y_true = idct(f(A0).*dct(b));

	runs = 1;	% multiple runs to get average times
	for run = 1:runs
		disp(run);
		para.spacing = 1;
		% para.tol = 1e-6;
		% para.truesol = y_true;
		para.method = 'standard';
		fprintf("%s:\t", para.method);
		tic;
		[y_std, hist_std, info_std] = arnoldi_matfun(A, b, f, m, para);
		times(testidx, 1) = times(testidx, 1) + toc/runs;
		fprintf("\n");
		nerr_std{testidx} = vecnorm(hist_std.sol - y_true)./norm(y_true);
		iter_std{testidx} = hist_std.iter;


		para.method = 'rand';
		fprintf("%s:\t\t", para.method);
		tic;
		[y_rand, hist_rand, info_rand] = arnoldi_matfun(A, b, f, m, para);
		times(testidx, 2) = times(testidx, 2) + toc/runs;
		fprintf("\n");
		nerr_rand{testidx} = vecnorm(hist_rand.sol - y_true)./norm(y_true);
		iter_rand{testidx} = hist_rand.iter;

		para.method = 'oblique';
		para.lastorth_method = 'chol';
		fprintf("%s:\t", para.method);
		tic;
		[y_obl, hist_obl, info_obl] = arnoldi_matfun(A, b, f, m, para);
		times(testidx, 3) = times(testidx, 3) + toc/runs;
		fprintf("\n");
		nerr_obl{testidx} = vecnorm(hist_obl.sol - y_true)./norm(y_true);
		iter_obl{testidx} = hist_obl.iter;

		para.method = 'oblique';
		para.lastorth_method = 'lsqr';
		para.lsqr_tol = 1e-12;
		fprintf("%s:\t", para.method);
		tic;
		[y_obl2, hist_obl2, info_obl2] = arnoldi_matfun(A, b, f, m, para);
		times(testidx, 4) = times(testidx, 4) + toc/runs;
		fprintf("\n");
		nerr_obl2{testidx} = vecnorm(hist_obl2.sol - y_true)./norm(y_true);
		iter_obl2{testidx} = hist_obl2.iter;

		para.method = 'oblique';
		para.lastorth_method = 'lsqr';
		para.lsqr_tol = 1e-1;
		fprintf("%s:\t", para.method);
		tic;
		[y_obl3, hist_obl3, info_obl3] = arnoldi_matfun(A, b, f, m, para);
		times(testidx, 5) = times(testidx, 5) + toc/runs;
		fprintf("\n");
		nerr_obl3{testidx} = vecnorm(hist_obl3.sol - y_true)./norm(y_true);
		iter_obl3{testidx} = hist_obl3.iter;

	end

end
	clear hist_std hist_rand hist_obl hist_obl2 hist_obl3;
	clear info_std info_rand info_obl info_obl2 info_obl3;
	clear y_std y_rand y_obl y_obl2 y_obl3;
	save("data_conv");


%% Figures:

for jj = 1:numtests

	leg = strings(0);
	figure;
	semilogy(iter_std{jj}, nerr_std{jj}, 'k-', 'LineWidth',2);	hold on;
	leg(end+1) = "Arnoldi";

	semilogy(iter_rand{jj}, nerr_rand{jj}, 'b-', 'LineWidth',2);	hold on;
	leg(end+1) = "randomized Arnoldi";

	semilogy(iter_obl{jj}, nerr_obl{jj}, 'r-', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (chol)";

	semilogy(iter_obl2{jj}, nerr_obl2{jj}, '--', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (lsqr $10^{-12}$)";
	
	semilogy(iter_obl3{jj}, nerr_obl3{jj}, 'm-.', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (lsqr $10^{-1}$)";

	legend(leg, 'fontsize', 22, 'location', 'northeast', 'interpreter', 'latex','box','off');
	xlim([-20, 320]);
	yy = ylim;
	yy(2) = 1e2;
	ylim(yy);
	xlabel('iteration','FontSize',20);
	ylabel('relative error','FontSize',20);
	titlestring = sprintf("Error - %s", testlabels(jj));
	title(titlestring,'FontSize',26);

	drawnow;

	ax = gca;
	ax.XAxis.FontSize = 20;
	ax.YAxis.FontSize = 20;

	filename = sprintf('fig/matfun_conv_%s.eps', testlabels(jj));
	exportgraphics(gcf, filename, 'ContentType', 'vector');

	% plot error ratio with respect to arnoldi:
	leg = strings(0);
	figure;

	semilogy(iter_rand{jj}, nerr_rand{jj}./nerr_std{jj}, 'b-', 'LineWidth',2);	hold on;
	leg(end+1) = "randomized Arnoldi";

	semilogy(iter_obl{jj}, nerr_obl{jj}./nerr_std{jj}, 'r-', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (chol)";

	p3 = semilogy(iter_obl2{jj}, nerr_obl2{jj}./nerr_std{jj}, '--', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (lsqr $10^{-12}$)";
	
	p4 = semilogy(iter_obl3{jj}, nerr_obl3{jj}./nerr_std{jj}, '-.', 'LineWidth',1.5);	hold on;
	leg(end+1) = "SRR-Arnoldi (lsqr $10^{-1}$)";
	p3.Color = p4.Color;
	p4.Color = 'm';

	legend(leg, 'fontsize', 22, 'location', 'northwest', 'interpreter', 'latex','box','off');
	xlim([-20, 320]);
	xlabel('iteration','FontSize',20);
	ylim([0.5, 6e2])
	ylabel('error ratio','FontSize',20);
	titlestring = sprintf("Error ratio - %s", testlabels(jj));
	title(titlestring,'FontSize', 26);
	
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
	