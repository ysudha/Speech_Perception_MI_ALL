function clus = spectral_clusters(mat, ch_list)
nch = size(mat,1);

ch = 1:nch;

D = diag(sum(mat));
L = D-mat; % definition of laplacian
[eigvec, eigvalues] = eig(L,D);

%k = 3; %number of clusters
eigvec = eigvec(:,1:k); %first k eigenvectors
idx = kmeans(eigvec,k); %running k means, where data points are the rows of the eigen vector matrix
%indices idx determine the cluster to which that data point belongs

b_ordered=[];
% clusters=[];
for i = 1:k
    %find cluster number i
    clus_k = find(idx==i);
    b_ordered= [b_ordered; clus_k];
%     clusters = [clusters; i*ones(length(clus_k),1)];
end

figure('Position', [0, 0, 1000, 1000]);
set(gca,'fontsize',6);
chnames = data.ch_names(b_ordered,:);
xlabels = cellstr(chnames);
imagesc(1:nch,1:nch,b(b_ordered,b_ordered))
set(gca,'XTick',[1:nch],'XTickLabel',xlabels);
xticklabel_rotate
set(gca,'YTick',[1:nch],'YTickLabel',xlabels);
title(str)
%print(get(get(gca, 'Title'),'String'),'-djpeg');
