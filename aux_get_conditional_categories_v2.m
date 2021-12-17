function [Lnew] = aux_get_conditional_categories_v2(L)

%posture_categories = {'other'};
posture_categories = {'prone','crawl_posture','standing','supine','sitting','other'};
posture_names = L(1).posture_names;
movement_categories = {'macro_still','proto_movement','elementary_movement','proficient_movement'};
%movement_categories = {'elementary_movement','proficient_movement'};
movement_names = L(1).movement_names;

Lnew = L;

new_movement_names = {};
ii = 1;
for i = 1:length(movement_names)
    if ~sum(strcmp(movement_names{i}, movement_categories))
        new_movement_names{ii} = movement_names{i};
        ii = ii+1;
    else
        for j = 1:length(posture_categories)
            %str1 = [movement_names{i} '-' posture_categories{j}];
            str1 = [posture_categories{j} '-' movement_names{i}];
            new_movement_names{ii} = str1;
            ii = ii+1;
        end
    end
end
% if length(posture_categories) == 4
% ORDER = [1,9,13,17, 2,10, 14, 18, 3, 11, 15, 19, 4, 12, 16, 20, 21, 5, 6, 7, 8];
% elseif length(posture_categories) == 5
%     ORDER = [1,10,15,20, 2,11,16,21, 3,12,17,22, 4,13,18,23, 5,14,19,24,25,6,7,8,9];
% 
% end
ORDER = 1:length(new_movement_names);
new_movement_names = new_movement_names(ORDER);
Ncats = length(new_movement_names);


for iRec = 1:length(L)
    %mask = logical(1-L(iRec).mask);
    mask = logical(ones(size(L(iRec).posture,1)));
    %[~, mov] = max(L(iRec).movement_oh,[],2);
    mov = L(iRec).movement_predicted+1;
    %mov_oh = L(iRec).movement;
    
    mov_named = movement_names(mov);
    
    %[~, pos] = max(L(iRec).posture,[],2);
    pos = L(iRec).posture_predicted+1;
    pos_oh = cat2onehot(pos,length(posture_names));
    pos_named = posture_names(pos); 
    
    % Convert movements to new classes
    N = length(mask);
    new_mov_named = convert_movement_classes(pos_named, mov_named, movement_categories, posture_categories);

    
    for i = 1:N
        mov(i) = find(strcmp(new_mov_named{i},new_movement_names));
    end
         
    mov_oh = cat2onehot(mov,Ncats);
    [mov_oh, new_movement_names2] = fuse_movement_categories(mov_oh, new_movement_names);
    Lnew(iRec).movement = mov_oh;
    Lnew(iRec).movement_names = new_movement_names2;
    try
    [pos_oh, new_posture_names] = fuse_posture_categories(pos_oh, posture_names);
    catch
        disp('');
    end
    Lnew(iRec).posture = pos_oh;
    Lnew(iRec).posture_names = new_posture_names;
 
    Lnew(iRec).movement_distribution = sum(Lnew(iRec).movement,1);
    Lnew(iRec).movement_distribution_norm = Lnew(iRec).movement_distribution/sum(Lnew(iRec).movement_distribution);
    
    Lnew(iRec).posture_distribution = sum(Lnew(iRec).posture,1);
    Lnew(iRec).posture_distribution_norm = Lnew(iRec).posture_distribution/sum(Lnew(iRec).posture_distribution);   
    
    
%     Imov = [1,4];
%     Lnew(iRec).movement_distribution(:,Imov) = [];
%     Lnew(iRec).movement_distribution(:,Imov) = [];
%     Lnew(iRec).movement_names(Imov) = [];
end



end


function [mov_oh, movement_names] = fuse_movement_categories(mov_oh, movement_names)
    I1 = find(strcmp(movement_names,'pivot_L'));
    I2 = find(strcmp(movement_names,'pivot_R'));
    I3 = find(strcmp(movement_names,'turn_L'));
    I4 = find(strcmp(movement_names,'turn_R'));
    
    mov_oh(:,I1) = mov_oh(:,I1) + mov_oh(:,I2);
    mov_oh(:,I3) = mov_oh(:,I3) + mov_oh(:,I4);
    
    mov_oh(:,[I2, I4]) = [];
    
    movement_names{I1} = 'pivot';
    movement_names{I3} = 'turn';
    movement_names([I2,I4]) = [];
end

function [pos_oh, posture_names] = fuse_posture_categories(pos_oh, posture_names)
    I1 = find(strcmp(posture_names,'side_L'));
    I2 = find(strcmp(posture_names,'side_R'));
    pos_oh(:,I1) = pos_oh(:,I1) + pos_oh(:,I2);
    pos_oh(:,I2) = [];
    posture_names{I1} = 'side';
    posture_names(I2) = [];
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
                new_pos_str = 'other';
            end
            %mov_named{i} = [mov_str '-' new_pos_str];
            mov_named{i} = [new_pos_str '-' mov_str];
              
        else
            mov_named{i} = mov_str;

        end

    end
    
    out = mov_named;

end



function out = cat2onehot(R, Ncats)
    out = zeros(length(R(:)),Ncats);
    
    for i = 1:length(R)
        out(i,R(i)) = 1;
    end
end