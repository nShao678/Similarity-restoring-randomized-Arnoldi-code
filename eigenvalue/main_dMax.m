addpath('../utils/');
clear all
close all
rng(1)
n0 = 10000;
n00 = 4;
numeval = 10;
n = n0*n00+numeval;
A0 = zeros(n,1);
for ii = 1:n00
    A0((ii-1)*n0+1:ii*n0) = 10^(ii)+10^(ii-1)*randn(n0,1);
end
A0(n0*n00+1:n0*n00+numeval) = randn(numeval,1);


dMinSet = [numeval,1.5*numeval,2*numeval,2.5*numeval];
dMaxSet = [3*numeval,4*numeval,5*numeval,6*numeval];

time = zeros(length(dMinSet),length(dMaxSet),3,2);
hist = cell(length(dMinSet),length(dMaxSet),3,2);

flag = zeros(length(dMinSet),length(dMaxSet),3,2);
for iterdMin = 1:length(dMinSet)	% size
    for iterdMax = 1:length(dMaxSet)
        fprintf("%d ", iterdMax);
        for jj = 1:2	% testcase
            if mod(jj,2) == 0
                A1 = randn(n-1,1);	% nonsym
            else
                A1 = zeros(n-1,1);	% symmetric
            end
            A = @(x) fun(x,A0,A1);
            q = randn(n,1);
            para.hist = 1;
            para.matvecMax = 10000;
            d = 100;
            Omega = sparse_sign(d,n,8);
            para.dMin = dMinSet(iterdMin);
            para.dMax = dMaxSet(iterdMax);
            para.sketch.Omega = @(x) Omega*x;
            para.sketch.d = d;
            para.type = 'smallestreal';
           
            para.getinfo = false;	% get additional info, but slower runtime
            para.method = 'KS';
            para.ks_reorth = true;	% use reorth in standard KS
            tic
            [u,matvec,hist{iterdMin,iterdMax,1,jj},~,flag(iterdMin,iterdMax,1,jj)] = KS(A,q,numeval,para);
            time(iterdMin,iterdMax,1,jj) = toc;

            para.method = 'OKS';
            para.lastorth_method = 'chol';	% 'chol' or 'lsqr'
            tic
            [ocu,ocmatvec,hist{iterdMin,iterdMax,2,jj},~,flag(iterdMin,iterdMax,2,jj)] = KS(A,q,numeval,para);
            time(iterdMin,iterdMax,2,jj) = toc;
%
            
            para.method = 'RKS';
            tic
            [ru,rmatvec,hist{iterdMin,iterdMax,3,jj},~,flag(iterdMin,iterdMax,3,jj)] = KS(A,q,numeval,para);
            time(iterdMin,iterdMax,3,jj) = toc;

        end
        
    end
    fprintf("\n");
end
clear A0 A1 Omega
clear ocu q ru u 
save('datadMax')
%% heatmap for time
figure
cmin = min(time(:));
cmax = max(time(:));
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
 
for jj = 1:2
    for ii = 1:3
        nexttile;
        imagesc(time(:,:,ii,jj));
        switch ii
            case 1
                title('KS');
            case 2
                title('SRR-KS')
            case 3
                title('RKS')
        end
        for iii = 1:length(dMinSet)
            for jjj = 1:length(dMaxSet)
                text(jjj, iii, num2str(floor(time(iii,jjj,ii,jj))), ...
            'HorizontalAlignment','center', ...
            'Color','w','FontWeight','bold');
                if flag(iii,jjj,ii,jj)==0
                    w = 0.4;   
                    h = 0.3;  
                    line([jjj-w, jjj+w], [iii-h, iii-h], 'Color', 'w', 'LineWidth', 1); % bottom
                    line([jjj-w, jjj+w], [iii+h, iii+h], 'Color', 'w', 'LineWidth', 1); % top
                    line([jjj-w, jjj-w], [iii-h, iii+h], 'Color', 'w', 'LineWidth', 1); % left
                    line([jjj+w, jjj+w], [iii-h, iii+h], 'Color', 'w', 'LineWidth', 1); % right
                end
            end
        end
        xticklabels(string(dMinSet));
        yticklabels(string(dMaxSet));
        axis equal tight;
        clim([cmin cmax])
    end
end

annotation('textbox', [0.15 0.95 0.7 0.05], ...
    'String', 'Hermitian matrix', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 14);
annotation('line', [0.06 0.9], [0.53 0.53], 'Color', 'k', 'LineWidth', 1.5);
annotation('textbox', [0.15 0.47 0.7 0.05], ...
    'String', 'Non-Hermitian matrix', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 14);
colormap jet;
cb = colorbar;
cb.Layout.Tile = 'east';
set(gcf, 'Color', 'w');


%% history for bad case

figure
hold on
ii = 1;
jj = 3;
res = vecnorm(hist{ii,jj,1,1}.res);
plot(res,'k-o','LineWidth',2,'DisplayName','KS')
res = vecnorm(hist{ii,jj,2,1}.res);
plot(res,'r-','LineWidth',2,'DisplayName','SRR-KS')
res = vecnorm(hist{ii,jj,3,1}.res);
plot(res,'b-d','LineWidth',2,'DisplayName','RKS')
hold off
set(gca,'yscale','log')
set(gcf, 'Color', 'w');
legend('fontsize',18,'Location','southwest','box','off')

figure
hold on
ii = 1;
jj = 3;
res = vecnorm(hist{ii,jj,1,2}.res);
plot(res,'k-o','LineWidth',2,'DisplayName','KS')
res = vecnorm(hist{ii,jj,2,2}.res);
plot(res,'r-','LineWidth',2,'DisplayName','SRR-KS')
res = vecnorm(hist{ii,jj,3,2}.res);
plot(res,'b-d','LineWidth',2,'DisplayName','RKS')
hold off
set(gca,'yscale','log')
set(gcf, 'Color', 'w');
legend('fontsize',18,'Location','southeast','box','off')
%%

function y = fun(x,A0,A1)
y = fft(x);
y = A0.*y+[A1.*y(1:end-1);0];
y = ifft(y);
end

