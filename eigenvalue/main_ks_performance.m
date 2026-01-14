clear all
addpath('../utils/');

rng(1)
n0 = 500;
time = zeros(5,8,4);
matvec = time;
for ii = 1:5	% size
    ii
    for jj = 1:8	% testcase
        fprintf("%d ", jj);
        n = n0*2^ii;
        a = linspace(2,10,n)';
        switch ceil(jj/2)
            case 1
                A0 = exp(a/10);
            case 2
                A0 = log(a+1);
            case 3
                A0 = 1+1./(a.^2);
            case 4
                A0 = 0.99.^a;
        end
        if mod(jj,2) == 0
            A1 = randn(n-1,1);	% nonsym
        else
            A1 = zeros(n-1,1);	% symmetric
        end
        A = @(x) fun(x,A0,A1);
        q = randn(n,1);
        para.hist = 1;
        para.matvecMax = 20000;
        numeval = 20;
        d = 8*numeval;
        Omega = sparse_sign(d,n,8);
        para.sketch.Omega = @(x) Omega*x;
        para.sketch.d = d;

        para.getinfo = 0;	% get additional info, but slower runtime


        para.method = 'KS';
        para.ks_reorth = true;	% use reorth in standard KS
        tic
        [u,matvec(ii,jj,1),hist] = KS(A,q,numeval,para);
        time(ii,jj,1) = toc;

        para.method = 'OKS';
        para.lastorth_method = 'chol';	% 'chol' or 'lsqr'
        tic
        [ocu,matvec(ii,jj,2),ochist] = KS(A,q,numeval,para);
        time(ii,jj,2) = toc;

        para.lastorth_method = 'lsqr';	% 'chol' or 'lsqr'
        tic
        [olu,matvec(ii,jj,3),olhist] = KS(A,q,numeval,para);
        time(ii,jj,3) = toc;

        para.method = 'RKS';
        tic
        [ru,matvec(ii,jj,4),rhist] = KS(A,q,numeval,para);
        time(ii,jj,4) = toc;

    end
    fprintf("\n");
end
clear A0 A1 Omega a q
clear u ocu olu ru
save('dataperformance')
%%
xx = n0*2.^(1:5);


for jj = 1:8


    figure
    hold on
    plot(xx, time(:,jj,1),'k-o','LineWidth',2,'DisplayName','KS')
    plot(xx, time(:,jj,2),'r-','LineWidth',2,'DisplayName','SRR-KS (chol)')
    plot(xx, time(:,jj,3),'g-^','LineWidth',2,'DisplayName','SRR-KS (lsqr)')
    plot(xx, time(:,jj,4),'b-d','LineWidth',2,'DisplayName','RKS')
    hold off
    ax = gca;
    xLimits = [min(xx),max(xx)]; % [xmin xmax]

    x_pos = 10^(log10(xLimits(2)) - 0.15*(log10(xLimits(2)) - log10(xLimits(1))));
    y_pos = mean(time(1,jj,:));

    % Add the text
    text(x_pos, y_pos, ['$f_{',num2str(ceil(jj/2)),'}$'], ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', ...
        'Color','k', ...
        'FontWeight','bold', ...
        'FontSize',20, ...
        'Interpreter','latex');
    legend('fontsize',18,'Location','northwest','box','off')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gcf, 'Color', 'w');

end


function y = fun(x,A0,A1)
y = fft(x);
y = A0.*y+[A1.*y(1:end-1);0];
y = ifft(y);
end

