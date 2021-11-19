% Low - BS >0 and High - BS>0 edges are good. knn subgraphs of energy


chlist = [1:25 28];%[1 3 4 8 12 13 2 6 7 11 15 16];%[1:25 28];
 %% (1) load graphs, edge = highVOT
chlist = [1:25 28];
miVOT789 = load('mi_matrix_VOT789_smallwin.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 13;
nch = length(chlist);
ntim = [[1:25:51] [101:25:326]];
ntim = ntim*2-202;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT
diff_highVOT = zeros(np,nch,nch,nwin-3);

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
    for iwin = 4:nwin
       baseline = (mi789{1,1} + mi789{2,1})./2;
       mat2 = squeeze(mi789{iwin,1} - baseline); 
       mat = squeeze(mat2);
       diff_highVOT(ipat,:,:,iwin-3) = mat(chlist,chlist);
    end

end

nwin = nwin-3;
ntim = ntim(4:end);
diff_highVOT(diff_highVOT<0)=0; %< is the right one
diff_highVOT = abs(diff_highVOT);

high_VOT_energy = zeros(np,nwin);

for iwin = 1:nwin
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
% load('ChanNamesandLocs.mat');
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
load('ChanNamesandLocs.mat');
chnames = chanNamesLocs(:,1);
chnames = chnames(chlist);
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
    for iwin = 1:nwin
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
    
inode = 13;
iwin = 1;
ipat = np;
figure;
for iwin = 1:nwin
mat = squeeze(diff_highVOT(ipat,:,:,iwin));
mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});

chl = nearestNeighbors{inode};

empty_graph = zeros(length(chlist),length(chlist));
empty_graph(chl,chl) = mat;

G = graph(empty_graph,chnames);
%G = graph(empty_graph);%,chnames);
plot(G,'XData',cell2mat(chanNamesLocs(chlist,2)),'YData',cell2mat(chanNamesLocs(chlist,3)))
az=-90;
ele = 90;
view([az ele])
waitforbuttonpress;
end


