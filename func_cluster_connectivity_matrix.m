function[b_ordered,box_borders]= func_cluster_connectivity_matrix(con_matrix, str,ch_names,k)

% this function taken in a connectivity matrix, using MI or conditional mi,
% a str for the title of the figure and a vector of channel names in the
% original matrix.
%it outputs a clustered figure.

% returns - indices of the clustered chanels, with the end of each cluster
% given in box borders
% 
data.ch_names = ch_names;

b = con_matrix;

nch = size(b,1);

%% clustering
%list of vertices/channels
ch = 1:nch;

bD = diag(sum(b'));
bL = bD-b;

[U, d] = eig(k, m, 'vector');
[d, ind] = sort(d);
U = U(:, ind);

[eigvec eigvalues] = eig(bL,bD);

%k = 3; %number of clusters
eigvec = eigvec(:,1:k); %first k eigenvectors
idx = kmeans(eigvec,k); %running k means, where data points are the rows of the eigen vector matrix
%indices idx determine the cluster to which that data point belongs

b_ordered=[];
% clusters=[];
box_borders = zeros(1,k);
for i = 1:k
    %find cluster number i
    clus_k = find(idx==i);
    b_ordered= [b_ordered; clus_k];
    box_borders(i) = b_ordered(end);
end

for i = 1:k
   box_borders(i) =  find(b_ordered==box_borders(i));
end

figure('Position', [0, 0, 1000, 1000]);
set(gca,'fontsize',6);
chnames = data.ch_names(b_ordered,:);
xlabels = cellstr(chnames);
imagesc(1:nch,1:nch,b(b_ordered,b_ordered))
set(gca,'XTick',[1:nch],'XTickLabel',xlabels);
%xticklabel_rotate
set(gca,'YTick',[1:nch],'YTickLabel',xlabels);
title(str)
%print(get(get(gca, 'Title'),'String'),'-djpeg');

end
