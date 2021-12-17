function figure_4B()

load('example_data.mat');
S = aux_get_conditional_categories_v2(S);

% GET DISTRIBUTIONS
P = zeros(length(S),length(S(1).posture_distribution_norm));
M  = zeros(length(S),length(S(1).movement_distribution_norm));

N_babies = size(P,1);

for k = 1:N_babies
    P(k,:) = S(k).posture_distribution_norm;
    M(k,:) = S(k).movement_distribution_norm;
end

D = [P, M];

IDs = {'1','2','3'};
% AGE CLASSIFIER
ages = [S(:).age]';
ages(ages > 16.01) = 16.01;
age_estimate = BIMS_LOSO(ages,D,[],1);

cc = [0.725, 0.608, 0.82];
figure();

aux_plot_regression2(ages, age_estimate,cc,'Age classifier','Pearson','linear',IDs);

xticks([4, 8, 12, 16]); xticklabels({'4','8','12','16+'});
yticks([4, 8, 12, 16]); yticklabels({'4','8','12','16+'});
xlabel('True age (months)');
ylabel('Predicted age (months)')
yl = get(gca,'ylim');
yyaxis right;
ylim(yl);
yticks([4,8,12,16]); yticklabels({'0','33','66','100'});
ylabel('BIMS');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
title('BIMS classifier')

% AIMS CLASSIFIER
aims = [S(:).aims_score]';
aims_estimate = BIMS_LOSO(aims,D,[],1);

figure();
aux_plot_regression2(aims, aims_estimate,cc,['AIMS Classifier'],'Pearson','linear',IDs);
xlabel('True AIMS score (total)');
ylabel('Predicted AIMS score (total)')
title('AIMS classifier')

end




function x_estimate =  BIMS_LOSO(xvals, D, Iexc, delta_grid)

% Get posture and movement distributions

if max(xvals) > 30 % AIMS case
    x_grid = min(xvals):delta_grid:max(xvals);
else
    x_grid = 10:delta_grid:16; % 4:1:16 used in study
end

x_estimate = zeros(size(D,1),1);
N_babies = size(D,1);
% LOSO loop
for subject_to_test = 1:N_babies
    
    % Define set of subjects in training set
    subjects_to_train = setxor(subject_to_test,1:N_babies);
    
    % Number of features to model
    featdim = size(D,2);
    
    % mean and variance matrices (feature x age group)
    mu = zeros(featdim,length(x_grid));
    sigma2 = zeros(featdim,length(x_grid));
    
    x = 1;
    for xval = x_grid
        
        % Age distance of each baby from the current age grid point
        d = abs(xval-xvals);
        
        % Maximum allowed age difference
        maxdist = 1*delta_grid;
        
        % Find babies within range, and exclude test subject
        x_i = find(d <= maxdist);
        x_i = intersect(x_i,subjects_to_train);
        Nmin = 1; % 3 used in study
        if length(x_i) < Nmin
            [~, I] = sort(d,'ascend');
            x_i = I;
            x_i(x_i == subject_to_test) = [];
            %disp('n');
            x_i = x_i(1:Nmin);
        end
        
        
        % Calculate means and variances
        for feat = 1:featdim
            ii = D(x_i,feat);
            m = mean(ii);
            mu(feat,x) = m;
            sigma2(feat,x) = max([std(ii),0.0001]);
        end
        
        x = x+1;
    end
    
    
    % Take held-out test subject and test each feature against each age
    % grid point
    
    likelihood = zeros(featdim,length(x_grid));
    for xval = 1:length(x_grid)
        for feat = 1:featdim
            likelihood(feat,xval) = normpdf(D(subject_to_test,feat),mu(feat,xval),(sigma2(feat,xval)));
        end
    end
    likelihood = likelihood + 0.000001;
    post_log = exp(nansum(log(likelihood),1));
    
    post_log = post_log./sum(post_log);
    
    x_estimate(subject_to_test) = sum(post_log.*x_grid);
    
end


end



function aux_plot_regression2(X, Y, color1, title_str, ctype,fitType,IDs)
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

Mdl = fitlm(X,Y,fitType);
minn = min([floor(min(X)), floor(min(Y))]);
maxx = max([ceil(max(X)), ceil(max(Y))]);
ddd = 0.5;
xvec = linspace(minn-ddd, maxx+ddd, 100)';
[ypred1, yint1] = predict(Mdl, xvec);
plot(xvec,ypred1,'Color',color1,'LineWidth',2); hold on;
plot3(X,Y,1:length(X),'x','Color',color1);

pstr = ['(p = ' num2str(p,1) ')'];





% Get unique IDs & Plot multiple recordings
uqID = unique(IDs);
for i = 1:length(uqID)
    Is = find(strcmp(IDs,uqID{i}));
    if length(Is) > 1
        plot3(X(Is),Y(Is),Is,'.--','Color',1-(1-color1).^0.75,'LineWidth',1.5,'MarkerSize',25);
    end
    
    
end


patch([xvec;flipud(xvec)],[yint1(:,1);flipud(yint1(:,2))],color1,'FaceAlpha',0.2,'EdgeColor','none');
title_str(1) = upper(title_str(1));
title_str = strrep(title_str,'_',' ');


xlim([minn-ddd,maxx+ddd]);
ylim([minn-ddd,maxx+ddd]);

xx = xlim;
yy = ylim;
cval = round(c,2);
if cval == 1
    cstr = '0.999';
else
    cstr = num2str(cval);
end
text(xx(1)+ 0.05*(xx(2)-xx(1)), yy(2)- 0.01*(yy(2)-yy(1)),['r = ' cstr ' ' pstr],'FontSize',12,'VerticalAlignment','top');


end
