function figure_4A()
%close all;
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

% Age-group models
ages = [S(:).age];
[muP, sigmaP, x_grid] = gaussianModel(ages,P,1);
[muM, sigmaM, ~] = gaussianModel(ages,M,1);

% Rearrange movement names & order:
movement_names = {'turn','transition','prone-proto','pivot','prone-elementary','prone-proficient','crawl_posture-proto','crawl_posture-elementary','crawl_posture-proficient','standing-proto','standing-elementary','standing-proficient'};

for i = 1:length(movement_names)
    try
        II(i) = find(contains(S(1).movement_names,movement_names{i}));
    catch
        disp('err')
    end
end
muM = muM(:,II); sigmaM = sigmaM(:,II); M = M(:,II);
movement_names = {'Roll','Transition','Prone-Proto','Pivot','Prone-Elementary','Prone-Fluent','Crawl-Proto','Crawl-Elementary','Crawl-Fluent','Standing-Proto','Standing-Elementary','Standing-Fluent'};
%movement_names = {'Roll','Transition','Proto','Pivot','Elementary','Fluent','Proto','Elementary','Fluent','Proto','Elementary','Fluent'};


% Plot posture distributions
POSTURE_ORDER = [2,1,3,4,5,6];
posture_names = {'Supine','Prone','Side','Crawl','Sitting','Standing'};
figure();
plot(nan,nan);hold on;
cm = colormap('lines'); cm = cm(1:6,:); cm = cm([1,2,6,4,5,3],:);
plotViolin(x_grid,muP(:,POSTURE_ORDER),sigmaP(:,POSTURE_ORDER), ages, P(:,POSTURE_ORDER), posture_names,0,cm,1:6);
xlabel('Age (months)')
xpos = 4.5; yp=6+0.05;
plot([xpos,xpos],[yp,yp+0.25],'k');
plot([xpos,xpos],[yp,yp+0.25],'k_');
text(xpos+0.3,yp+0.125,'25%','FontSize',14);
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Plot movement distributions
figure();
plot(nan,nan);hold on;
movement_groups = {[1,2], [3,4,5,6], [7,8,9], [10,11,12]};
mm = [2,4,3,3];
dd = 0.5;
ypos = [];
ii = 0;
for i = 1:length(mm)
    for j = 1:mm(i)
        if j == 1
            ypos = [ypos, ii+2*dd];
            ii = ypos(end);
        else
            ypos = [ypos, ii+dd];
            ii = ypos(end);
        end
    end
end

rr = 0.5;
cm2 = [cm(1,:); cm(1,:).^rr; cm(2,:); cm(2,:).^(rr); cm(2,:).^(rr^2); cm(2,:).^(rr^3); cm(4,:); cm(4,:).^(rr); cm(4,:).^(rr^2); cm(6,:); cm(6,:).^(rr); cm(6,:).^(rr^2)];
plotViolin(x_grid,muM,sigmaM, ages, M, movement_names,0,cm2,ypos);
xpos = 4.5; yp=ypos(end)+0.15;
plot([xpos,xpos],[yp,yp+0.25],'k');
plot([xpos,xpos],[yp,yp+0.25],'k_');
text(xpos+0.3,yp+0.125,'25%','FontSize',14);
xlabel('Age (months)')

end

function [mu, sd, range] = pdfNormFit(x,y)
%determine from cumulative distribution;
ycdf = cumsum(y);
val0 = 0.5;
val1 = 0.25; val2=0.75;
[vvdist, vvinds] = sort(abs(ycdf-val0),'ascend');
vvinds = sort(vvinds(1:2));
v = vvdist(vvinds);
vratio = abs(v(1))/sum(abs(v));
ages = [x(vvinds(1)), x(vvinds(2))];
mu = (1-vratio)*ages(1) + (vratio)*ages(2);


[~, i1] = min(abs(ycdf-val1));
[~, i2] = min(abs(ycdf-val2));
range = [x(i1), x(i2)];

sd = (abs(mu-x(i1)) + abs(mu-x(i2)))/2;

end

function M = compress(M)
tval = 0.2;
d = M-tval;
I = find(d > 0);
M(I) = tval + log10(d(I)+1)/3;

end

function patchplot(x, y,i, c)
scale = 1;
patch([x; flipud(x)], [scale*y+i; scale*flipud(-y)+i],c,'EdgeColor','none','FaceAlpha',0.95);
end

function distplot(mu, range, i)
d = 0.15;

plot([range(1), range(2)],[i+d,i+d],'k','Marker','.','MarkerSize',10,'LineWidth',1);
plot(mu,i+d,'xr');hold on;

end

function plotViolin(x_grid, mu, sigma, ages, yvals, names,do_scaling,cm,ypos)
N = length(names);
scale = (max(mu,[],1));
if do_scaling
    mu_scaled = compress(mu);
    maxval = max(mu(:));
    mu_scaled = mu_scaled./maxval*0.5;
    yvals = compress(yvals)/maxval * 0.5 ;
    
else
    sc = 0.5*0.85;
    maxval = max(mu(:));
    mu_scaled = mu./maxval * sc;
    yvals = yvals/maxval * sc;
end

mu_norm = mu./sum(mu,1);
for i = 1:size(mu,2)
    [mean_ages(i), std_ages(i), range_ages(i,:)] =pdfNormFit(x_grid, mu_norm(:,i));
end

for i = 1:N
    ii = ypos(i);
    patchplot(x_grid, mu_scaled(:,i),ii,cm(i,:)); hold on;
    plot(x_grid,ii*ones(size(x_grid)),'k');hold on;
    
end

for i = 1:N
    ii = ypos(i);
    distplot(mean_ages(i), range_ages(i,:),ii);
    % scatter(ages,yvals(:,i) + ii, 10,cm(i,:).^0.5); % Plot individual points
end

xlim([x_grid(1), x_grid(end)]);
xticks([4,8,12,16,20]);
yticks(ypos);
yticklabels(names);
set(gca,'TickLabelInterpreter','none');
end


function [mu, sigma2, x_grid] =  gaussianModel(xvals, D, delta_grid)


% Get posture and movement distributions
x_grid = 10:1:16; % 4:1:16 is used in the paper
N_babies = size(D,1);

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
    maxdist = delta_grid;
    
    % Find babies within range
    x_i = find(d <= maxdist);
    Nmin = 1; % 3 is used in paper
    if length(x_i) < Nmin
        [~, I] = sort(d,'ascend');
        x_i = I;
        x_i = x_i(1:Nmin);
    end
    
    
    
    % Calculate  mean and variance
    for feat = 1:featdim
        ii = D(x_i,feat);
        m = mean(ii);
        %m = median(ii);
        mu(feat,x) = m;
        sigma2(feat,x) = max([std(ii),0.0001]);
    end
    
    x = x+1;
end
mu =mu'; sigma2 = sigma2'; x_grid = x_grid';




end
