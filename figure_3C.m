function figure_3C()
close all;
% Load data
load('example_data.mat');
L = S;
labels_A = L(1).posture_names;
labels_B = L(1).movement_names;

% Fuse L/R categories
[LB_fused,~] = fuse_mov(L,labels_B);
[LA_fused,~] = fuse_pos(L,labels_A);

labels_A_fused = {'Prone','Supine','Side','Crawl','Sitting','Standing'};
labels_B_fused = {'Still','Roll','Pivot','Proto','Elementary','Fluent','Transition'};


% Regression analysis
regression_analysis(L,LA_fused,labels_A_fused);
regression_analysis(L,LB_fused,labels_B_fused);

% Bland-Altman analysis
ba_analysis(L,LA_fused,labels_A_fused);
ba_analysis(L,LB_fused,labels_B_fused);


end

function out = cat2onehot_s(vec,Ncats)
out = zeros(length(vec),Ncats);
for i = 1:length(vec)
    out(i,vec(i)) = 1;
end
end

function out = onehot2cat_s(M)
out = zeros(size(M,1),1);
for i = 1:size(M,1)
    [~, idx] = max(M(i,:));
    out(i) = idx;
end
end

function [L_fused, labels_fused] = fuse_mov(L, labels)
% Side L/R
to_fuse = [2,3; 4,5]; to_del = [3,5];
L_fused = L;
for i = 1:length(L)
    CM = L(i).movement_confmat;
    pp = cat2onehot_s(L(i).movement_predicted+1,length(labels));
    tt = cat2onehot_s(L(i).movement_target+1,length(labels));
    for j = 1:size(to_fuse,1)
        id1 = to_fuse(j,1);
        id2 = to_fuse(j,2);
        pp(:,id1) = pp(:,id1) + pp(:,id2);
        tt(:,id1) = tt(:,id1) + tt(:,id2);
        CM(id1,:) = CM(id1,:) + CM(id2,:);
        CM(:,id1) = CM(:,id1) + CM(:,id2);
    end
    pp(:,to_del) = []; tt(:,to_del) = [];
    CM(:,to_del) = []; CM(to_del,:) = [];
    L_fused(i).preds = onehot2cat_s(pp) - 1;
    L_fused(i).targets = onehot2cat_s(tt) - 1;
    L_fused(i).ConfMat = CM;
end
labels_fused = labels;
labels_fused(to_del) = [];
end

function [L_fused, labels_fused] = fuse_pos(L, labels)
% Side L/R
to_fuse = [3,4]; to_del = [4];
L_fused = L;
for i = 1:length(L)
    CM = L(i).posture_confmat;
    pp = cat2onehot_s(L(i).posture_predicted+1,length(labels));
    tt = cat2onehot_s(L(i).posture_target+1,length(labels));
    for j = 1:size(to_fuse,1)
        id1 = to_fuse(j,1);
        id2 = to_fuse(j,2);
        pp(:,id1) = pp(:,id1) + pp(:,id2);
        tt(:,id1) = tt(:,id1) + tt(:,id2);
        CM(id1,:) = CM(id1,:) + CM(id2,:);
        CM(:,id1) = CM(:,id1) + CM(:,id2);
    end
    pp(:,to_del) = []; tt(:,to_del) = [];
    CM(:,to_del) = []; CM(to_del,:) = [];
    L_fused(i).preds = onehot2cat_s(pp) - 1;
    L_fused(i).targets = onehot2cat_s(tt) - 1;
    L_fused(i).ConfMat = CM;
end
labels_fused = labels;
labels_fused(to_del) = [];
end

function regression_analysis(L, L2, category_names)

meansTar = zeros(length(L),length(category_names));
meansPred = zeros(length(L),length(category_names));
%figure();
for i = 1:length(L2)
    pred = L2(i).preds +1; mask = logical(1-L(i).mask);
    tar = L2(i).targets +1; pred = pred(mask); tar = tar(mask);
    Ntot = length(tar);
    for j = 1:length(category_names)
        Ntar = sum(tar == j);
        Npred = sum(pred == j);
        meansTar(i,j) = Ntar/Ntot;
        meansPred(i,j) = Npred/Ntot;
    end
end

figure();
for i = 1:length(category_names)
    subplot(3,ceil(length(category_names)/3),i);
    
    aux_plot_regression2(100*meansTar(:,i),100*meansPred(:,i), [0.725, 0.608, 0.82], category_names{i}, 'Pearson','Linear',1)
    
end

end



function out = ba_analysis(L, L2, category_names)

meansTar = zeros(length(L),length(category_names));
meansPred = zeros(length(L),length(category_names));
%figure();
for i = 1:length(L2)
    pred = L2(i).preds +1; mask = logical(1-L(i).mask);
    tar = L2(i).targets +1; pred = pred(mask); tar = tar(mask);
    Ntot = length(tar);
    for j = 1:length(category_names)
        Ntar = sum(tar == j);
        Npred = sum(pred == j);
        meansTar(i,j) = Ntar/Ntot;
        meansPred(i,j) = Npred/Ntot;
    end
end
ages = [L.age]';
figure();
for i = 1:length(category_names)
    subplot(3,ceil(length(category_names)/3),i);
    
    aux_plot_ba2(100*meansTar(:,i),100*meansPred(:,i), [0.725, 0.608, 0.82], category_names{i},ages)
    
end

end






function aux_plot_ba2(X, Y, color1, title_str, ages)
if nargin < 3
    color1 = 'b';
    title_str = '';
    ctype = 'Spearman';
    fitType = 'linear';
elseif nargin <4
    title_str = '';
    ctype = 'Spearman';
    fitType= 'linear';
elseif nargin <5
    ctype = 'Spearman';
    fitType = 'linear';
elseif nargin <6
    fitType = 'linear';
end
Yo = Y;
Y = Y-X;

I1 = ~isnan(X); I2 = ~isnan(Y); Itot = I1&I2;
X = X(Itot); Y = Y(Itot);

diffSTD = std(Y);
diffMean = mean(Y);

% t-test
[~,diffMeanP, ~, TSTATS] = ttest(Y,0);
RPC = 1.96*diffSTD;
% [kspH, KSP] = kstest((Y-diffMean)/diffSTD);
pstr = ['(p = ' num2str(diffMeanP,1) ')'];

% Plot difference data & mean
plot(X,Y,'x','Color',color1); hold on;
xvec = [min(X), max(X)]';
plot(xvec,[diffMean, diffMean],'Color',color1,'LineWidth',2);

% Fit regression model to distribution data
Mdl = fitlm(ages,Yo,'linear');

% Compute slope
[ypred1, yint1] = predict(Mdl, xvec);
vv = (ypred1(end)-ypred1(1))/(xvec(end)-xvec(1));
vv = abs(vv);
% Plot dashed lines
plot(xvec,[0+vv,0+vv],'--','Color',color1,'LineWidth',1.5);
plot(xvec,[0-vv,0-vv],'--','Color',color1,'LineWidth',1.5);
sensAxis = [0-vv, 0+vv];
pSens = mean(Y > sensAxis(1) & Y < sensAxis(2));

pstr = ['t(' num2str(TSTATS.df) ') = ' num2str(TSTATS.tstat,2) ' ' pstr];

% plot 95% confidence interval of data
yci = [diffMean + 1.96*diffSTD, diffMean + 1.96*diffSTD; diffMean - 1.96*diffSTD, diffMean - 1.96*diffSTD];
patch([xvec;flipud(xvec)],[yci(1,:)'; yci(2,:)'],color1,'FaceAlpha',0.2,'EdgeColor','none');

title_str(1) = upper(title_str(1));
title_str = strrep(title_str,'_',' ');

xlim([xvec(1),xvec(2)+eps]);
ylim([-10,10]);

xx = xlim;
yy = ylim;
text(xx(1)+ 0.05*(xx(2)-xx(1)), yy(2)- 0.01*(yy(2)-yy(1)),[pstr],'FontSize',14,'VerticalAlignment','top');

xticks([0, 25, 50, 75, 100]); %yticklabels({'0%', '25%', '50%','75%', '100%'});
try
    yticks([-10, round(yci(2,1),1), 0, round(yci(1,1),1), 10]); %xticklabels({'-10%','0%','10%'});
catch
    
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
title(title_str);

end

function aux_plot_regression2(X, Y, color1, title_str, ctype,fitType,plotAxis)
if nargin < 3
    color1 = 'b';
    title_str = '';
    ctype = 'Spearman';
    fitType = 'linear';
elseif nargin <4
    title_str = '';
    ctype = 'Spearman';
    fitType= 'linear';
elseif nargin <5
    ctype = 'Spearman';
    fitType = 'linear';
elseif nargin <6
    fitType = 'linear';
end

I1 = ~isnan(X); I2 = ~isnan(Y); Itot = I1&I2;
[c, p] = corr(X(I1&I2),Y(I1&I2),'type',ctype);
plot(X,Y,'x','Color',color1); hold on;
%scatter(X, Y,color1,'FaceAlpha',0.4); hold on;
xvec = linspace(floor(min(X)), ceil(max(X)), 100)';
Mdl = fitlm(X,Y,fitType);

%Mdl = fitlm(X,Y,'linear');
[ypred1, yint1] = predict(Mdl, xvec);

%if p < 0.001
%    pstr = '(p < 0.001)';
%else
pstr = ['(p = ' num2str(p,1) ')'];
%end

plot(xvec,ypred1,'Color',color1,'LineWidth',2);
%if strcmp(fitType,')
%c = sqrt(Mdl.Rsquared.Ordinary);
%p = Mdl.coefTest;
patch([xvec;flipud(xvec)],[yint1(:,1);flipud(yint1(:,2))],color1,'FaceAlpha',0.2,'EdgeColor','none');
%title([title_str 'R-squared=' num2str(c^2,2) ', p=' num2str(p,1) ', N=' num2str(sum(Itot))], 'interpreter', 'none');
title_str(1) = upper(title_str(1));
title_str = strrep(title_str,'_',' ');
title(title_str, 'interpreter','none')
minn = min([floor(min(X)), floor(min(Y))]);
maxx = max([ceil(max(X)), ceil(max(Y))]);
xlim([minn,maxx]);
ylim([minn,maxx]);

%xlim([0, 100]);
%ylim([0, 100]);
%ylim([ymin, ymax]);
%ylim([-2.5, 100]);
xx = xlim;
yy = ylim;
cval = round(c,2);
if cval == 1
    cstr = '0.999';
else
    cstr = num2str(cval);
end
text(xx(1)+ 0.05*(xx(2)-xx(1)), yy(2)- 0.01*(yy(2)-yy(1)),['r = ' cstr ' ' pstr],'FontSize',14,'VerticalAlignment','top');
%text(xx(1)+ 0.05*(xx(2)-xx(1)), yy(2)- 0.01*(yy(2)-yy(1)),['R^2=' num2str(c^2,2)],'FontSize',12,'VerticalAlignment','top');
%if 0%plotAxis
%xticks([0, 25, 50, 75, 100]); %yticklabels({'0%', '25%', '50%','75%', '100%'});
xticks([]);
yticks([0, 25, 50, 75, 100]); %xticklabels({'-10%','0%','10%'});
%else
%    yticks([]);
%end
%grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',14)
end
