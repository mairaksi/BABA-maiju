function figure_2C()
close all;
load('example_data.mat');
L = S;

catA = {'Prone','Supine','Side L','Side R','Crawl posture','Sitting','Standing'};
catB = {'Still','Roll L','Roll R','Pivot L','Pivot R','Proto','Elementary','Fluent','Transition'};


% ANNOTATION VS AGE:
ages = [L(:).age];
I = ~isnan(ages);
ages = ages(I)';

[Dpos, Dmov] = getDistributions(L(I));

plotCategoryCorrelations(Dpos, ages, catA, 'Age (months)', 'Proportion (%)',3);
sgtitle('Posture track age correlation')

plotCategoryCorrelations(Dmov, ages, catB, 'Age (months)', 'Proportion (%)',3);
sgtitle('Movement track age correlation')

% ANNOTATION VS AIMS:
aims = [L(:).aims_score];
I = ~isnan(aims);
aims = aims(I)';

[Dpos, Dmov] = getDistributions(L(I));
plotCategoryCorrelations(Dpos, aims, catA, 'AIMS (total)', 'Proportion (%)',3);
sgtitle('Posture track AIMS correlation')

plotCategoryCorrelations(Dmov, aims, catB, 'AIMS (total)', 'Proportion (%)',3);
sgtitle('Movement track AIMS correlation')




end


function out = plotCategoryCorrelations(D, xvar, names, xlab, ylab,Nrows,fname)
figure();
N = size(D,2);
if nargin < 6
    Nrows = 2;
end
Ncols = ceil(N/Nrows);

for i = 1:N
    subplot(Nrows,Ncols,i)
    
    % Settings for regression plot
    fitType = 'quadratic';
    base_color = [0.725, 0.608, 0.82]; % Base color for plots

    % Plot regression
    aux_plot_regression2(xvar,100*D(:,i), base_color, [names{i}],'Pearson',fitType);
    xlabel(xlab);
    if mod((i-1),Nrows) == 0
        ylabel(ylab);
    end
    
end



end

function aux_plot_regression2(X, Y, color1, title_str, ctype,fitType)
if nargin < 3
    color1 = 'b';
    title_str = '';
    ctype = 'Pearson';
    fitType = 'linear';
elseif nargin <4
    title_str = '';
    ctype = 'Pearson';
    fitType= 'linear';
elseif nargin <5
    ctype = 'Pearson';
    fitType = 'linear';
elseif nargin <6
    fitType = 'linear';
end

I1 = ~isnan(X); I2 = ~isnan(Y);

% Compute correlation
[c, p] = corr(X(I1&I2),Y(I1&I2),'type',ctype);
plot(X,Y,'x','Color',color1); hold on;
xvec = linspace(floor(min(X)), ceil(max(X)), 100)';
Mdl = fitlm(X,Y,fitType);
[ypred1, yint1] = predict(Mdl, xvec);

plot(xvec,ypred1,'Color',color1,'LineWidth',2);
patch([xvec;flipud(xvec)],[yint1(:,1);flipud(yint1(:,2))],color1,'FaceAlpha',0.2,'EdgeColor','none');

% Plot formatting
title_str(1) = upper(title_str(1));
title_str = strrep(title_str,'_',' ');
title(title_str, 'interpreter','none')
xlim([floor(min(X)), ceil(max(X))]);
ymax = ceil(max([max(Y), max(yint1(:,2))])); ymin = floor(min([min(Y), min(yint1(:,1))]));
yts = [0, 50, 100];
if ymin < -2.5
    ymin = -2.5;
end
dd = yts-ymax; I = dd>0;
if sum(I) == 0
    ymax = 100;
else
    yts = yts(I); dd=dd(I); 
    II = find(min(dd));
    ymax = yts(II);
end

ylim([ymin, ymax]);
ylim([-2.5, 100]);
xx = xlim;
yy = ylim;
pstr = ['(p = ' num2str(p,1) ')'];
text(xx(1)+ 0.05*(xx(2)-xx(1)), yy(2)- 0.01*(yy(2)-yy(1)),['r = ' num2str(c,2) ' ' pstr],'FontSize',12,'VerticalAlignment','top');

yticks([0, 50, 100]); yticklabels({'0%', '50%', '100%'});

end


function [Dpos, Dmov] = getDistributions(L)
% Compute posture and movement distributions from data
Dpos = zeros(length(L), size(L(1).posture_oh,2));
Dmov  = zeros(length(L), size(L(1).movement_oh,2));


for i = 1:length(L)
    mask = logical(1-L(i).mask);
    pos = L(i).posture_oh(mask,:);
    mov = L(i).movement_oh(mask,:);
    Dpos(i,:) = mean(pos,1);
    Dmov(i,:) = mean(mov,1);
end
end