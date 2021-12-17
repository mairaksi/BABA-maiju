function figure_3A()
load('example_data.mat');
L = S;

% Gather data
posture_names = L(1).posture_names;
movement_names = L(1).movement_names;
cpc_encs = []; mask = [];
for i = 1:length(L)
    mask = [mask; L(i).mask];
    cpc_encs = [cpc_encs; L(i).cpc_encs];
end
mask = logical(1-mask);

% Compute tsne
T = tsne(cpc_encs(mask,:));

% Plot tsne
pvec = [L.posture_target]' +1; pvec = pvec(mask);
mvec = [L.movement_target]' +1; mvec = mvec(mask);
[pvec, posture_names] = fuse_cats(pvec,posture_names,[3,4]);
[mvec, movement_names] = fuse_cats(mvec,movement_names,[2,3;4,5]);
[mvec, movement_names] = fuse_cats(mvec,movement_names,[7,2;5,3]);
figure();
aux_plot_tsne2(T,pvec,posture_names);

figure();
aux_plot_tsne2(T,mvec,movement_names);



end


function aux_plot_tsne2(Y, labs, category_names)
% Create custom colormap order
cm1 = colormap('lines');
cm1 = cm1(1:11,:);
cm1(8,:) = [0,0,0];
cm1(9,:) = [0.5, 0.5, 0.5];
cm1(10,:) = [1,1,0];
cm1(11,:) = [1,0,1];

% Maximum number of data points shown
N = 100000;
if N > size(Y,1)
    N = size(Y,1);
end

if sum(contains(category_names,'prone')) % Custom options for posture track style
    cat_names_str = flip({'Supine','Prone','Side','Crawl','Sitting','Standing'});
    ORD = flip([2,1,3,4,5,6]);
    category_names = category_names(ORD);
    labs = cat2onehot_s(labs,length(ORD));
    labs = labs(:,ORD);
    labs = onehot2cat_s(labs);
    
    cm1 = flipud(cm1([1,2,6,4,5,3],:));
elseif sum(contains(category_names,'macro_still')) % Custom options for Movement track style
    ORD = ([1,2,3,4,5]);
    category_names = category_names(ORD);
    labs = cat2onehot_s(labs,length(ORD));
    labs = labs(:,ORD);
    labs = onehot2cat_s(labs);
    cm1(5,:) = cm1(6,:);
    cm1 = cm1(ORD,:);
    cat_names_str = ({'Still','Proto','Elementary','Fluent','Transition'});
    cat_names_str = cat_names_str(ORD);
end

perm = randperm(size(Y,1));
perm = perm(1:N);

CC = cm1(labs,:);
hGS = scatter(Y(perm,1),Y(perm,2), 6, CC(perm,:), 'filled');
hold on;
for i = 1:length(category_names)
    h(i)=plot(nan,nan,'.','MarkerSize',20,'Color',cm1(i,:));
    
end
hh = legend(h,cat_names_str,'Location','SouthWest','interpreter','none','FontSize',14);
axis off; axis tight;

set(gca,'LooseInset',get(gca,'TightInset'));
legend boxoff;
pp = hh.Position;
set(hh,'Position',[0.05,0.05,pp(3),pp(4)]);
%xticks([]); yticks([]);

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

function [vec_fused, labels_fused] = fuse_cats(vec, labels, to_fuse)
% Side L/R
%to_fuse = [2,3; 4,5]; %to_del = [3,5];
to_del = to_fuse(:,2)';

vv = cat2onehot_s(vec,length(labels));
for j = 1:size(to_fuse,1)
    id1 = to_fuse(j,1);
    id2 = to_fuse(j,2);
    vv(:,id1) = vv(:,id1) + vv(:,id2);
end
vv(:,to_del) = [];
vec_fused = onehot2cat_s(vv);

labels_fused = labels;
labels_fused(to_del) = [];
end

