function figure_2B_irr()

close all;
clear;


load('example_data.mat');
L = S;

labels_A = {'Prone','Supine','Side L','Side R','Crawl','Sitting','Standing'};
labels_B = {'Still','Roll L','Roll R','Pivot L','Pivot R','Proto','Elementary','Fluent','Transition'};


%L = aux_get_conditional_categories_v2(L);
% Posture:
catA = labels_A;
catB = labels_B;

Ptot1 = []; Ptot2 = []; Mtot1 = []; Mtot2 = []; Mtot11 = []; Mtot22 = [];
Np = length(catA);
for i = 1:length(L)

    Nannot = size(L(i).posture,2);
    if Nannot > 1
        for j1 = 1:Nannot
            for j2 = 1:Nannot
                if j1 >= j2
                    continue
                end
                mask = logical(1-L(i).mask);
                pos1 = L(i).posture(:,j1);
                pos2 = L(i).posture(:,j2);
                mov1 = L(i).movement(:,j1);
                mov2 = L(i).movement(:,j2);
                
                mask2 = logical((mov1 > 0) & (mov2 > 0) & (pos1 > 0 & pos1 <= Np) & (pos2 > 0 & pos2 <= Np)) ;
                mask = mask & mask2;
                
                Ptot1 = [Ptot1; pos1(mask)];
                Ptot2 = [Ptot2; pos2(mask)];

                if sum(mov1(mask) == 0) > 0
                    disp('mov1');
                elseif sum(mov2(mask) == 0) > 0
                    disp('mov2');
                end
                [mov_cond1, cat_cond] = get_conditional_categories(pos1(mask),mov1(mask),catA,catB);
                mov_cond2 = get_conditional_categories(pos2(mask),mov2(mask),catA,catB);
                %cat_cond = labels_B;
                Mtot11 = [Mtot11; L(i).movement(mask,j1)];
                Mtot22 = [Mtot22; L(i).movement(mask,j2)];
                Mtot1 = [Mtot1; mov_cond1'];
                Mtot2 = [Mtot2; mov_cond2'];
            end
        end
    end
end

INDS_P = []; INDS_M = [];
for i = 1:length(Ptot1)
    if Ptot1(i) == 0 || Ptot2(i) == 0
        INDS_P = [INDS_P, i];
    end
    if Mtot1(i) == 0 || Mtot2(i) == 0
        INDS_M = [INDS_M, i];
    end
end
Ptot1(INDS_P) = [];
Ptot2(INDS_P) = [];
Mtot1(INDS_M) = [];
Mtot2(INDS_M) = [];
Mtot11(INDS_M) = [];
Mtot22(INDS_M) = [];


kappa_p = fleisskappa([Ptot1, Ptot2]);
kappa_m = fleisskappa([Mtot11, Mtot22]);
kappa_m2 = fleisskappa([Mtot1, Mtot2]);

disp(['Posture kappa=' num2str(kappa_p)]);
disp(['Movement kappa=' num2str(kappa_m)]);
disp(['Posture-conditional movement kappa=' num2str(kappa_m2)]);


% Get confusion matrices

CMa = zeros(length(catA));
CMb2 = zeros(length(catB ));
CMb = zeros(length(cat_cond));
for i = 1:length(Ptot1)
    val1 = Ptot1(i); val2=Ptot2(i);
    if val1 <= length(catA) && val2 <= length(catA)
      CMa(Ptot1(i),Ptot2(i)) = CMa(Ptot1(i),Ptot2(i)) + 1;
    end

end
for i = 1:length(Mtot1)
    %    val1 = Mtot1(i); val2=Mtot2(i);
    %if val1 <= length(catB) && val2 <= length(catB)
    CMb(Mtot1(i),Mtot2(i)) = CMb(Mtot1(i),Mtot2(i)) + 1;
    %end
end

for i = 1:length(Mtot11)
    %    val1 = Mtot1(i); val2=Mtot2(i);
    %if val1 <= length(catB) && val2 <= length(catB)
    CMb2(Mtot11(i),Mtot22(i)) = CMb(Mtot11(i),Mtot22(i)) + 1;
    %end
end
CMa = CMa + CMa';
CMb = CMb + CMb';
CMb2 = CMb2 + CMb2';

figure();
subplot(1,2,1)
plotConfMatrix(CMa,catA)
xlabel('Annotator 1');
ylabel('Annotator 2');
 xtickangle(30);
 ytickangle(30);
 H = get(gca,'Title');
tit = H.String;
 title(['Inter-rater agreement' newline 'Posture track' newline 'kappa=' num2str(kappa_p,3) ', ' tit])
 set(gca,'TickLabelInterpreter','none')
subplot(1,2,2)
plotConfMatrix(CMb2,catB)
xlabel('Annotator 1');
ylabel('Annotator 2');
 xtickangle(30);
 ytickangle(30);
 H = get(gca,'Title');
tit = H.String;
 title(['Inter-rater agreement' newline 'Movement track' newline 'kappa=' num2str(kappa_m,3) ', ' tit])
 set(gca,'TickLabelInterpreter','none')
figure();

plotConfMatrix(CMb,cat_cond);
xlabel('Annotator 1');
ylabel('Annotator 2');
 xtickangle(30);
 ytickangle(30);
 H = get(gca,'Title');
tit = H.String;
 title(['Inter-rater agreement' newline 'Posture conditional Movement track' newline 'kappa=' num2str(kappa_m2,3) ', ' tit])
 set(gca,'TickLabelInterpreter','none')

end



function plotConfMatrix(varargin)
switch (nargin)
    case 0
        confmat = 1;
        labels = {'1'};
    case 1
        confmat = varargin{1};
        labels = 1:size(confmat, 1);
        titl = '';
    case 2
        confmat = varargin{1};
        labels = varargin{2};
        titl = '';
    otherwise
        confmat = varargin{1};
        labels = varargin{2};
        titl = varargin{3};
end

for i = 1:length(labels)
    
    str = labels{i};
    str(1) = upper(str(1));
    str= strrep(str,'_',' ');
    labels{i} = str;
end

confmat(isnan(confmat))=0; % in case there are NaN elements
numlabels = size(confmat, 1); % number of labels
confmat = confmat';
% calculate the percentage accuracies
%confpercent = 100*confmat./repmat(sum(confmat, 1),numlabels,1);
confpercent = 100*confmat./repmat(sum(confmat, 2),1,numlabels);
confpre = 100*confmat./repmat(sum(confmat, 1),numlabels,1);
confrec = 100*confmat./repmat(sum(confmat, 2),1,numlabels);
confpercent2 = 2.* (confpre .* confrec) ./ (confpre + confrec);
cp2=diag(confpercent2);
for i=1:numlabels
    disp([labels{i} ' F1: ' num2str(cp2(i),2)]);
end
disp(['UWF1:' num2str(mean(cp2))]);
% plotting the colors
imagesc(confpercent);
title(sprintf([titl ' Accuracy: %.1f%%, UWF1: %.1f%%'], 100*trace(confmat)/sum(confmat(:)), mean(cp2(~isnan(cp2)))));
%title(sprintf([titl ' Accuracy: %.2f%%'], 100*trace(confmat)/sum(confmat(:))),'FontSize',16);

%xlabel('Predicted class'); ylabel('Target class');

% set the colormap
%colormap(flipud(gray));
colormap(flipud(bone));

% Create strings from the matrix values and remove spaces
textStrings = num2str([confpercent(:), confmat(:)], '%.0f%%\n%d\n');
%textStrings = num2str([confpercent(:)], '%.0f%%');
textStrings = strtrim(cellstr(textStrings));

% Create x and y coordinates for the strings and plot them
[x,y] = meshgrid(1:numlabels);
hStrings = text(x(:),y(:),textStrings(:), ...
    'HorizontalAlignment','center');

% Get the middle value of the color range
midValue = mean(get(gca,'CLim'));

% Choose white or black for the text color of the strings so
% they can be easily seen over the background color
textColors = repmat(confpercent(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2),'FontSize',8);

% Setting the axis labels
set(gca,'XTick',1:numlabels,...
    'XTickLabel',labels,...
    'YTick',1:numlabels,...
    'YTickLabel',labels,...
    'TickLength',[0 0]);

ytickangle(30);
xtickangle(30);
end



function [Bnew, new_movement_names] = get_conditional_categories(A, B, posture_names, movement_names)

posture_categories = {'Prone','Supine','Crawl','Sitting','Standing','Other'};
movement_categories = {'Still','Proto','Elementary','Fluent'};


new_movement_names = {};
ii = 1;

for i = 1:length(posture_categories)
    for j = 1:length(movement_names)
        if ~sum(strcmp(movement_names{j}, movement_categories))
            if ~any(strcmp(new_movement_names,movement_names{j}))
                new_movement_names{ii} = movement_names{j};
                ii = ii+1;
            end
        else
            str1 = [posture_categories{i} '-' movement_names{j}];
            new_movement_names{ii} = str1;
            ii = ii+1;
        end
    end
end
Ncats = length(new_movement_names);
B_named = movement_names(B);
A_named = posture_names(A);


% Convert movements to new classes
N = length(A_named);
new_B_named = convert_movement_classes(A_named, B_named, movement_categories, posture_categories);

for i = 1:N
    Bnew(i) = find(strcmp(new_B_named{i},new_movement_names));
end

end


function out = convert_movement_classes(pos_named, mov_named, movement_categories, posture_categories)
N = length(pos_named);
for i = 1:N
    pos_str = pos_named{i};
    mov_str = mov_named{i};
    
    if sum(strcmp(mov_str,movement_categories))
        if sum(strcmp(pos_str,posture_categories))
            new_pos_str = pos_str;
        else
            new_pos_str = 'Other';
        end
        mov_named{i} = [new_pos_str '-' mov_str];
        
    else
        mov_named{i} = mov_str;
        
    end
    
end

out = mov_named;

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

function kappa = fleisskappa(X)
% function kappa = fleisskappa(X)
%
% Computes inter-annotator agreement level for data in X, where each row of
% X corresponds to a sample and each column of X correspond to categorical
% decisions (X(i,j) = [1, 2,..,k] from each annotator.

N = size(X,1);
k = max(X(:));
n = size(X,2);

% Proportion of assignments to 0 and 1
pj = zeros(k,1);
for j = 1:k
    pj(j) = sum(X(:) == j)./(N*n);
end

nij = zeros(N,2);
for i = 1:N
    for j = 1:k
        nij(i,j) = sum(X(i,:) == j);
    end
end

coef = 1/(n*(n-1));
Pi = zeros(N,1);
for i = 1:N
    for j = 1:k        
        Pi(i) = Pi(i)+nij(i,j)*(nij(i,j)-1);
    end
    Pi(i) = Pi(i)*coef;
end

Pdash = 1/N*sum(Pi);

Pe = sum(pj.^2);

kappa = (Pdash-Pe)/(1-Pe);
end
