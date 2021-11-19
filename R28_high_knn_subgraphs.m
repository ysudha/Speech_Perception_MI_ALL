% Low - BS >0 and High - BS>0 edges are good. knn subgraphs of energy


chlist = [1:25 28];%[1 3 4 8 12 13 2 6 7 11 15 16];%[1:25 28];
 %% (1) load graphs, edge = highVOT
chlist = [1:25 28];
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
ntim = [101:5:301]*2-202;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT
diff_highVOT = zeros(np,nch,nch,nwin-2);

% remove the negative values after bias correction
for ipat = 1:np
    for iwin = 1:nwin
        mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1};
        mi789(mi789<0)=0;
        miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1} = mi789;
    end
end

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    for iwin = 3:nwin
       mat2 = squeeze(mi789{iwin,1} - mi789{1,1}); 
       mat = squeeze(mat2);
       diff_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

diff_highVOT(diff_highVOT<0)=0; %< is the right one
diff_highVOT = abs(diff_highVOT);

high_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_highVOT(ipat,:,:,iwin));
     %   bar(ipat,norm(mat,'fro'));
        high_VOT_energy(ipat,iwin) = norm(mat,'fro');
    end
end

%let's plot the energy of baseline normalized graphs
figure;
for i = 1:np-1
    plot(ntim,high_VOT_energy(i,:),'b-*','LineWidth',4);
    hold on;
end
plot(ntim,high_VOT_energy(np,:),'r-*','LineWidth',4);
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: (HighVOT- BS) >0');
%title('Energy of matrices: (LowVOT- BS) <0');

%% statistics for high
X1 = high_VOT_energy;

stats = zeros(1,nwin);
stats_h_bs = zeros(1,np);
stats_high_gt_bs = zeros(np,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('high_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_h_bs(ipat) = sum(stats<0.05);
    stats_high_gt_bs(ipat,:) = stats<0.05;
end

figure;
imagesc(ntim, 1:22,stats_high_gt_bs);
title('Significant values High-BS>0');
xlabel('Time');
ylabel('Subjects');
ax = gca;
ax.XGrid = 'on'
grid minor;
axis square;
colorbar;

%% knn subgraphs
% 
load('ChanNamesandLocs.mat');
% chnames = chanNamesLocs(:,1);
% chnames = chnames(chlist);
% 
% % G = graph(mask);%,chnames);
% % plot(G,'XData',cell2mat(chanNamesLocs(chlist,2)),'YData',cell2mat(chanNamesLocs(chlist,3)),'ZData',cell2mat(chanNamesLocs(chlist,4)))
% 
% imp_ch_centers = zeros(length(4:25),length(chlist));
% for k_n = 4:25
% %X = cell2mat(chanNamesLocs(1,2:4));
% X = cell2mat(chanNamesLocs([1:25 28],2:4));
% %figure; imagesc(squareform(pdist(X)));
% %colormap jet; colorbar;
% dMat = squareform(pdist(X));
% nearestNeighbors = cell(1,length(chlist));
% for inode = 1:length(chlist)
%     vec = dMat(inode,:);
%     [val,ind] = sort(vec);
%     nearestNeighbors{inode} = ind(1:k_n);
% end
%     
% 
% % find energy for each subplot
% 
% for inode = 1:length(chlist)
%     
%     low_VOT_energy = zeros(np,nwin);
%     for iwin = 1:41
%        % figure;hold on;
%         for ipat = 1:np
% 
%             mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
%             mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
%             low_VOT_energy(ipat,iwin) = norm(mat,'fro');
%         end
%     end
%     
% %     figure;
% %     subplot(2,1,1);
% %     for i = 1:np-1
% %         plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
% %         hold on;
% %     end
% %     plot(ntim,low_VOT_energy(np,:),'r-*','LineWidth',4);
% %     ylabel('Frobenius norm');
% %     xlabel('Time');
% %     str = sprintf('Energy of subgraph centered around %d',inode);
% %     title(str);
%     
%     %subplot(2,1,2);
%         X1 = low_VOT_energy;
% 
%     stats = zeros(1,nwin);
%     stats_l_bs = zeros(1,np);
%     stats_low_gt_bs = zeros(np,nwin);
% 
%     for ipat = 1:np
%         pat_list = 1:np;
%         pat_list(ipat)=[];
%         for iwin = 1:nwin
%             [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
%             stats(iwin) = p;
% 
%         end
%         stats = fdr(stats);
%         %fprintf('low_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
%         stats_l_bs(ipat) = sum(stats<0.05);
%         stats_low_gt_bs(ipat,:) = stats<0.05;
%     end
% 
% %     imagesc(ntim, 1:22,stats_low_gt_bs);
% %     str = sprintf('Significant values %d',sum(stats_low_gt_bs(end,:)));
% %     title(str);
% %     xlabel('Time');
% %     ylabel('Subjects');
% %     ax = gca;
% %     ax.XGrid = 'on'
% %     grid minor;
% %     axis square;
% %     colorbar;
% 
%     imp_ch_centers(k_n-3,inode)=sum(stats_low_gt_bs(end,:));
% end
% 
% 
% end
% 
% figure; imagesc(imp_ch_centers);
% 
% thresh_imp = imp_ch_centers;
% thresh_imp(thresh_imp<35)=0;
% 
% figure; imagesc(thresh_imp)

%% look at individual knn subgraphs

k_n = 4;
%X = cell2mat(chanNamesLocs(1,2:4));
X = cell2mat(chanNamesLocs([1:25 28],2:4));
%figure; imagesc(squareform(pdist(X)));
%colormap jet; colorbar;
dMat = squareform(pdist(X));
nearestNeighbors = cell(1,length(chlist));
for inode = 1:length(chlist)
    vec = dMat(inode,:);
    [val,ind] = sort(vec);
    nearestNeighbors{inode} = ind(1:k_n);
end
    

% find energy for each subplot

sig_windows_aphasia = zeros(length(chlist),nwin);
for inode = 1:length(chlist)
    
    high_VOT_energy = zeros(np,nwin);
    for iwin = 1:41
        for ipat = 1:np
            mat = squeeze(diff_highVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            high_VOT_energy(ipat,iwin) = norm(mat,'fro');
        end
    end
    
    figure;
    subplot(2,1,1);
    for i = 1:np-1
        plot(ntim,high_VOT_energy(i,:),'b-*','LineWidth',4);
        hold on;
    end
    plot(ntim,high_VOT_energy(np,:),'r-*','LineWidth',4);
    ylabel('Frobenius norm');
    xlabel('Time');
    str = sprintf('Energy of subgraph centered around %d',inode);
    title(str);
    
    subplot(2,1,2);
    X1 = high_VOT_energy;

    stats = zeros(1,nwin);
    stats_h_bs = zeros(1,np);
    stats_high_gt_bs = zeros(np,nwin);

    for ipat = 1:np
        pat_list = 1:np;
        pat_list(ipat)=[];
        for iwin = 1:nwin
            [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
            stats(iwin) = p;

        end
        stats = fdr(stats);
        %fprintf('low_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
        stats_h_bs(ipat) = sum(stats<0.05);
        stats_high_gt_bs(ipat,:) = stats<0.05;
    end

    imagesc(ntim, 1:22,stats_high_gt_bs);
    sig_windows_aphasia(inode,:)=stats_high_gt_bs(end,:);
    str = sprintf('Significant values %d',sum(stats_high_gt_bs(end,:)));
    title(str);
    xlabel('Time');
    ylabel('Subjects');
    ax = gca;
    ax.XGrid = 'on'
    grid minor;
    axis square;
    colorbar;

end

figure; %imagesc(sig_windows_aphasia)
[r, c] = size(sig_windows_aphasia);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_aphasia);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
%% fixed knn value and node, whats the subgraph plot?

k_n = 4;

nearestNeighbors = cell(1,length(chlist));
for inode = 1:length(chlist)
    vec = dMat(inode,:);
    [val,ind] = sort(vec);
    nearestNeighbors{inode} = ind(1:k_n);
end
    
inode = 18;
iwin = 1;
ipat = np;
figure;
for iwin = 1%:41
mat = squeeze(diff_highVOT(ipat,:,:,iwin));
mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});

chl = nearestNeighbors{inode};

empty_graph = zeros(length(chlist),length(chlist));
empty_graph(chl,chl) = mat;

G = graph(empty_graph,chnames);
%G = graph(empty_graph);%,chnames);
plot(G,'XData',cell2mat(chanNamesLocs(chlist,2)),'YData',cell2mat(chanNamesLocs(chlist,3)),'ZData',cell2mat(chanNamesLocs(chlist,4)))
az=-90;
ele = 90;
view([az ele])
waitforbuttonpress;
end

%% spectral clustering
%func_cluster_connectivity_matrix(con_matrix, str,ch_names,k)
healthy_matrices = squeeze(sum(diff_highVOT(1:np-1,:,:,:)));
aphasic_matrices = squeeze((diff_highVOT(np,:,:,:)));
new_edge_matrices = aphasic_matrices - healthy_matrices;
new_edge_matrices(new_edge_matrices<0)=0;


M.matrices = new_edge_matrices;

mat = (squeeze(M.matrices(:,:,1)));
n = numel(mat(:,1));
gamma = 0.8;
omega = 1;
T=length(squeeze(M.matrices(1,1,:)));
B=spalloc(n*T,n*T,n*n*T+2*n*T);
twomu=0;
for s=1:T
    k=sum(squeeze(M.matrices(:,:,s)));
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:n]+(s-1)*n;
    B(indx,indx)=squeeze(M.matrices(:,:,s))-gamma*k'*k/twom;
end
twomu=twomu+2*omega*n*(T-1);
B = B + omega*spdiags(ones(n*T,2),[-n,n],n*T,n*T);
[S,Q] = genlouvain(B);
Q = Q/twomu;
partition_multi = reshape(S,n,T);

% What does the partition of nodes into communities look like?

figure; imagesc(partition_multi)

ch = partition_multi(:,1);

%% lets plot the energies for these individual partitions

nbrainRegions = 4;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = find(ch==1);
nearestNeighbors{2} = find(ch==2);
nearestNeighbors{3} = find(ch==3);
nearestNeighbors{4} = find(ch==4);
nearestNeighbors{5} = find(ch==5);


sig_windows_aphasia = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
    
    low_VOT_energy = zeros(np,nwin);
    for iwin = 1:41
        for ipat = 1:np
            mat = squeeze(diff_highVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            low_VOT_energy(ipat,iwin) = norm(mat,'fro');
        end
    end
    figure;
    subplot(2,1,1);
    for i = 1:np-1
        plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
        hold on;
    end
    plot(ntim,low_VOT_energy(np,:),'r-*','LineWidth',4);
    ylabel('Frobenius norm');
    xlabel('Time');
    str = sprintf('Energy of subgraph centered around %d',inode);
    title(str);
    
    subplot(2,1,2);

    X1 = low_VOT_energy;

    stats = zeros(1,nwin);
    stats_l_bs = zeros(1,np);
    stats_low_gt_bs = zeros(np,nwin);

    for ipat = 1:np
        pat_list = 1:np;
        pat_list(ipat)=[];
        for iwin = 1:nwin
            [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
            stats(iwin) = p;

        end
        stats = fdr(stats);
        %fprintf('low_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
        stats_l_bs(ipat) = sum(stats<0.05);
        stats_low_gt_bs(ipat,:) = stats<0.05;
    end

    sig_windows_aphasia(inode,:)=stats_low_gt_bs(end,:);
 
    imagesc(ntim, 1:22,stats_low_gt_bs);
    str = sprintf('Significant values %d',sum(stats_low_gt_bs(end,:)));
    title(str);
    xlabel('Time');
    ylabel('Subjects');
    ax = gca;
    ax.XGrid = 'on'
    grid minor;
    axis square;
    colorbar;
    
end

figure; %imagesc(sig_windows_aphasia)
[r, c] = size(sig_windows_aphasia);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_aphasia);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

     %% plot this in the graphs

load('ChanNamesandLocs.mat');
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);
chl = nearestNeighbors{4};

empty_graph = zeros(length(chlist),length(chlist));
empty_graph1 = zeros(length(chl),length(chl));
%empty_graph(chl,chl) = mat;

G = graph(empty_graph,chnames);
G1 = graph(empty_graph1,chnamesOrig(chl));
%G = graph(empty_graph);%,chnames);
plot(G,'XData',cell2mat(chanNamesLocs(chlist,2)),'YData',cell2mat(chanNamesLocs(chlist,3)),'NodeFontSize',14); hold on;
p=plot(G1,'XData',cell2mat(chanNamesLocs(chl,2)),'YData',cell2mat(chanNamesLocs(chl,3)),'NodeFontSize',14)
p.NodeColor='r';
az=-90;
ele = 90;
view([az ele])
box off;

