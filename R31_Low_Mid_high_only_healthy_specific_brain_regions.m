% back to drawing board. How do the individual edges b,aseline normalized perform? Any
% differentiability?


chlist = [1:25 28];%[1 3 4 8 12 13 2 6 7 11 15 16];%[1:25 28];
 %% (1) load graphs, edge = lowVOT
chlist = [1:25 28];
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
ntim = [101:5:301]*2-202;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT
diff_lowVOT = zeros(np,nch,nch,nwin-2);

% remove the negative values after bias correction
for ipat = 1:np
    for iwin = 1:nwin
        mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1};
        mi123(mi123<0)=0;
        miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1} = mi123;
    end
end

for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    for iwin = 3:nwin
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       mat = squeeze(mat2);
       diff_lowVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

diff_lowVOT(diff_lowVOT<0)=0; %< is the right one
diff_lowVOT = abs(diff_lowVOT);

low_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
     %   bar(ipat,norm(mat,'fro'));
        low_VOT_energy(ipat,iwin) = norm(mat,'fro');
    end
end

%let's plot the energy of baseline normalized graphs
figure;
for i = 1:np-1
    plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
    hold on;
end
plot(ntim,low_VOT_energy(np,:),'r-*','LineWidth',4);
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: (LowVOT- BS) >0');
%title('Energy of matrices: (LowVOT- BS) <0');

%% statistics for low
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
    fprintf('low_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_l_bs(ipat) = sum(stats<0.05);
    stats_low_gt_bs(ipat,:) = stats<0.05;
end

figure;
imagesc(ntim, 1:22,stats_low_gt_bs);
title('Significant values Low-BS>0');
xlabel('Time');
ylabel('Subjects');
ax = gca;
ax.XGrid = 'on'
grid minor;
axis square;
colorbar;

 %% (2) load graphs, edge = highVOT
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
 %% (3)load graphs, edge = midVOT
chlist = [1:25 28];
miVOT456 = load('mi_matrix_VOT456_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
ntim = [101:5:301]*2-202;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT
diff_midVOT = zeros(np,nch,nch,nwin-2);
% remove the negative values after bias correction
for ipat = 1:np
    for iwin = 1:nwin
        mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1};
        mi456(mi456<0)=0;
        miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1} = mi456;
    end
end

for ipat = 1:np
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1}; %contains 43 time windows
    for iwin = 3:nwin
       mat2 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat = squeeze(mat2);
       diff_midVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

diff_midVOT(diff_midVOT<0)=0; %< is the right one
diff_midVOT = abs(diff_midVOT);
mid_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_midVOT(ipat,:,:,iwin));
     %   bar(ipat,norm(mat,'fro'));
        mid_VOT_energy(ipat,iwin) = norm(mat,'fro');
    end
end

%let's plot the energy of baseline normalized graphs
figure;
for i = 1:np-1
    plot(ntim,mid_VOT_energy(i,:),'b-*','LineWidth',4);
    hold on;
end
plot(ntim,mid_VOT_energy(np,:),'r-*','LineWidth',4);
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: (MidVOT- BS) >0');
%title('Energy of matrices: (LowVOT- BS) <0');

%% statistics for mid
X1 = mid_VOT_energy;

stats = zeros(1,nwin);
stats_m_bs = zeros(1,np);
stats_mid_gt_bs = zeros(np,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('mid_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_m_bs(ipat) = sum(stats<0.05);
    stats_mid_gt_bs(ipat,:) = stats<0.05;
end

figure;
imagesc(ntim, 1:22,stats_mid_gt_bs);
title('Significant values Mid-BS>0');
xlabel('Time');
ylabel('Subjects');
ax = gca;
ax.XGrid = 'on'
grid minor;
axis square;
colorbar;


%% let's look at individual indiffernces between low, mid and high for normals
ntim = ([101:5:301]-101)*2;

healthy_low_VOT_energy = low_VOT_energy(1:end-1,:);
healthy_mid_VOT_energy = mid_VOT_energy(1:end-1,:);
healthy_high_VOT_energy = high_VOT_energy(1:end-1,:);

figure; hold on;
plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
jbfill(ntim,mean(healthy_low_VOT_energy)+std(healthy_low_VOT_energy)*2.086/sqrt(21),mean(healthy_low_VOT_energy)-std(healthy_low_VOT_energy)*2.086/sqrt(21),'r','r',0,0.3);
hold on;
plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
jbfill(ntim,mean(healthy_high_VOT_energy)+std(healthy_high_VOT_energy)*2.086/sqrt(21),mean(healthy_high_VOT_energy)-std(healthy_high_VOT_energy)*2.086/sqrt(21),'b','b',0,0.3);
hold on;
plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
jbfill(ntim,mean(healthy_mid_VOT_energy)+std(healthy_mid_VOT_energy)*2.086/sqrt(21),mean(healthy_mid_VOT_energy)-std(healthy_mid_VOT_energy)*2.086/sqrt(21),'k','k',0,0.3);
hold on;
legend('Low','','High','','Mid','');

figure;
hold on;
ntim = ([101:5:301]-101)*2;
plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
legend('Low','High','Mid');
title('Energy of VOT condition (mean of all healthy people at each time)');
ylabel('Mean energy of MI matrices');
xlabel('Time (ms)');

%% comparing mid, low and high VOTs, average of healthy and aphasia guy
aphasia_low_VOT_energy = low_VOT_energy(end,:);
aphasia_mid_VOT_energy = mid_VOT_energy(end,:);
aphasia_high_VOT_energy = high_VOT_energy(end,:);

healthy_low_VOT_energy = low_VOT_energy(1:end-1,:);
healthy_mid_VOT_energy = mid_VOT_energy(1:end-1,:);
healthy_high_VOT_energy = high_VOT_energy(1:end-1,:);

figure; hold on;
plot(ntim,(aphasia_low_VOT_energy),'-.r','LineWidth',6);hold on;
plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
jbfill(ntim,mean(healthy_low_VOT_energy)+std(healthy_low_VOT_energy)*2.086/sqrt(21),mean(healthy_low_VOT_energy)-std(healthy_low_VOT_energy)*2.086/sqrt(21),'r','r',0,0.3);
hold on;
plot(ntim,(aphasia_high_VOT_energy),'-.b','LineWidth',6);hold on;
plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
jbfill(ntim,mean(healthy_high_VOT_energy)+std(healthy_high_VOT_energy)*2.086/sqrt(21),mean(healthy_high_VOT_energy)-std(healthy_high_VOT_energy)*2.086/sqrt(21),'b','b',0,0.3);
hold on;
plot(ntim,(aphasia_mid_VOT_energy),'-.k','LineWidth',6);hold on;
plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
jbfill(ntim,mean(healthy_mid_VOT_energy)+std(healthy_mid_VOT_energy)*2.086/sqrt(21),mean(healthy_mid_VOT_energy)-std(healthy_mid_VOT_energy)*2.086/sqrt(21),'k','k',0,0.3);
hold on;
legend('Low','','High','','Mid','');

figure;
hold on;
ntim = ([101:5:301]-101)*2;
plot(ntim,(aphasia_low_VOT_energy),'-.r','LineWidth',6);hold on;
plot(ntim,(aphasia_high_VOT_energy),'-.b','LineWidth',6);hold on;
plot(ntim,(aphasia_mid_VOT_energy),'-.k','LineWidth',6);hold on;

plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
legend('Low','High','Mid');
title('Energy of VOT condition (mean of all healthy people at each time) vs aphasic (dotted line)');
ylabel('Mean energy of MI matrices');
xlabel('Time (ms)');

%% comparing statistically contrasts at each time point

low_VOT = healthy_low_VOT_energy;
high_VOT = healthy_high_VOT_energy;
mid_VOT = healthy_mid_VOT_energy;
pval_mid_gt_low_and_high = zeros(1,nwin);
contrasts = zeros(np-1,nwin);
t_statistics = zeros(1,nwin);
for i = 1:nwin
    anovaMatrix = [low_VOT(:,i) mid_VOT(:,i)  high_VOT(:,i)];
   bb = normalize_multiple_patients(anovaMatrix);
   %bb = anovaMatrix;
    col = bb(:,2) - 0.5*(bb(:,1)+bb(:,3));
    contrasts(:,i) = col;
    [h,p,ci,stats] = ttest(col,0,'Tail','right');
    pval_mid_gt_low_and_high(i)=p;
    t_statistics(i) = stats.tstat;
end
figure;
plot(contrasts');hold on;
figure;
plot(mean(contrasts),'k');
plot(ntim,mean(contrasts),'-*k','LineWidth',4);hold on;
jbfill(ntim,mean(contrasts)+std(contrasts)*2.086/sqrt(21),mean(contrasts)-std(contrasts)*2.086/sqrt(21),'k','k',0,0.3);
hold on;
line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
xlabel('time'); ylabel('contrast vector');
figure;
bar(ntim,pval_mid_gt_low_and_high); hold on;
line([0,ntim(end)],[0.05 0.05],'Color','black','LineStyle','--');
ylabel('p-values');xlabel('starting point of time-windows')

figure;
bar(ntim,t_statistics); hold on;
line([0,ntim(end)],[1.725 1.725],'Color','black','LineStyle','--');
ylabel('t-statistic');xlabel('starting point of time-windows')
%pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

% %% comparing statistically at each time point, using ranova
% 
% low_VOT = healthy_low_VOT_energy;
% high_VOT = healthy_high_VOT_energy;
% mid_VOT = healthy_mid_VOT_energy;
% pval_mid_gt_low_and_high = zeros(1,nwin);
% for i = 1:nwin
%     anovaMatrix = [low_VOT(:,i) mid_VOT(:,i)  high_VOT(:,i)];
%    bb = normalize_multiple_patients(anovaMatrix);
% 
%    [p,tbl,stats] = anova2(bb);
%    %results = multcompare(stats)
%   
% 
%     pval_mid_gt_low_and_high(i)= p(1);
%     %fcdf(1./rr(1,4),rr(2,2),rr(1,2),'upper');
% end
% 
% figure;
% bar(ntim,pval_mid_gt_low_and_high); hold on;
% line([0,ntim(end)],[0.05 0.05],'Color','black','LineStyle','--');
% 
% %pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

max_t_cluster = sum(t_statistics(9:10));
%% let's see if we can do cluster analysis
[pval,thresh] = func_cluster_analysis(100,nwin,np,low_VOT,mid_VOT,high_VOT,max_t_cluster)
% 
% nperm = 500;
% max_tstats = zeros(1,nperm);
% for iperm = 1:nperm
%     pval_mid_gt_low_and_high = zeros(1,nwin);
%     contrasts = zeros(np-1,nwin);
%     t_statistics = zeros(1,nwin);
%     for i = 1:nwin
%         anovaMatrix = [low_VOT(:,i) mid_VOT(:,i)  high_VOT(:,i)];
%         bb = normalize_multiple_patients(anovaMatrix);
%         for j = 1:np-1
%             bb(j,:) = bb(j,randperm(3));
%         end
%         %bb = anovaMatrix;
%         col = bb(:,2) - 0.5*(bb(:,1)+bb(:,3));
%         contrasts(:,i) = col;
%         [h,p,ci,stats] = ttest(col,0,'Tail','right');
%         pval_mid_gt_low_and_high(i)=p;
%         t_statistics(i) = stats.tstat;
%     end
% %     figure;
% %     plot(contrasts');hold on;
% %     figure;
% %     %plot(mean(contrasts),'k');
% %     plot(ntim,mean(contrasts),'-*k','LineWidth',4);hold on;
% %     jbfill(ntim,mean(contrasts)+std(contrasts)*2.086/sqrt(21),mean(contrasts)-std(contrasts)*2.086/sqrt(21),'k','k',0,0.3);
% %     hold on;
% %     line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
% %     xlabel('time'); ylabel('contrast vector');
% %     figure;
% %     bar(ntim,pval_mid_gt_low_and_high); hold on;
% %     line([0,ntim(end)],[0.05 0.05],'Color','black','LineStyle','--');
% %     ylabel('p-values');xlabel('starting point of time-windows')
% %     
% %     figure;
% %     bar(ntim,t_statistics); hold on;
% %     line([0,ntim(end)],[1.76 1.76],'Color','black','LineStyle','--');
% %     ylabel('t-statistic');xlabel('starting point of time-windows')
% tmask = t_statistics>1.725;
% cluster_starts = strfind((t_statistics>1.725),[0 1])+1; %+1 as we dont want to include 0 in the string 011
% t_stats_clusters = zeros(1,length(cluster_starts));
% for j = 1:length(cluster_starts)
%     
%     within_cluster_incementer = cluster_starts(j);
%     while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
%         t_stats_clusters(j) = t_stats_clusters(j)+ t_statistics(within_cluster_incementer);
%         within_cluster_incementer = within_cluster_incementer+1;
%     end
% end
%   
% if(length(t_stats_clusters)>0)
%     max_tstats(iperm) = max(t_stats_clusters);
% else
%     max_tstats(iperm) = 0;
% end
% 
% end
% 
% figure;
% (hist(max_tstats));
% tstats_sorted = sort(max_tstats);
% thresh = tstats_sorted(round(nperm*0.95))
% 
% pval = sum(max_tstats>max_t_cluster)./nperm

%% let's look at knn subgraphs, for significance between mid, low and high 

load('ChanNamesandLocs.mat');

k_n = 26;
%X = cell2mat(chanNamesLocs(1,2:4));
X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 15;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = [3 4 8 9 13];
nearestNeighbors{2} = [12 14 17 18 23 ];
nearestNeighbors{3} = [10 11 15 16 20];
nearestNeighbors{4} = [21 22 26];
nearestNeighbors{5} = [1 2 5 6 7];
nearestNeighbors{6} = [19 24 25 26]; 
nearestNeighbors{7} = [8 9 13 18 17 22];

 nearestNeighbors{8} = [2     3     5     9    11    12    17    18    23];
 nearestNeighbors{9} =     [  6    20    21    22    24    25];
 nearestNeighbors{10} =     [13    14    16    19];
  nearestNeighbors{11} =     [     1     4     7     8    10    15    26];
    
  nearestNeighbors{12} = [1 2 3 4 5 6 7 9 10];
 nearestNeighbors{13} =     [ 11 15 16 19 20];
 nearestNeighbors{14} =     [8 13 17 18 12 14];
  nearestNeighbors{15} =     [21 22 23 24 25 26];
    
    
  
% find energy for each subplot, lowVOT, highVOT and midVOT

sig_windows_all = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
    
    high_VOT_energy_knn = zeros(np-1,nwin);
    low_VOT_energy_knn = zeros(np-1,nwin);
    mid_VOT_energy_knn = zeros(np-1,nwin);

    for iwin = 1:41
        for ipat = 1:np-1
            mat = squeeze(diff_highVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            high_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            low_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_midVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            mid_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
        end
    end
    
    pval_mid_gt_low_and_high = zeros(1,nwin);
    for i = 1:nwin
        anovaMatrix = [low_VOT_energy_knn(:,i) mid_VOT_energy_knn(:,i)  high_VOT_energy_knn(:,i)];
        bb = normalize_multiple_patients(anovaMatrix);
        col = 2*bb(:,2) - (bb(:,1)+bb(:,3));
        [h,p] = ttest(col,0,'Tail','right');
        pval_mid_gt_low_and_high(i)=p;
    end
   % pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

    sig_windows_all(inode,:) = pval_mid_gt_low_and_high;

end

sig_windows_all(sig_windows_all>0.05)=1;

figure; %imagesc(sig_windows_aphasia)
[r, c] = size(sig_windows_all);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_all);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
     

% %% new algorithm, delete a node at a time aphasia vs rest - lowVOT use max-sum tstat
% 
% load('ChanNamesandLocs.mat');
% 
% k_n = 26;
% %X = cell2mat(chanNamesLocs(1,2:4));
% X = cell2mat(chanNamesLocs([1:25 28],2:4));
%     
% % find energy for each subplot, lowVOT, highVOT and midVOT
% current_ch_list = 1:26;
% tstats_for_iterations = zeros(1,length(chlist));
% 
% % full netowork stats: no nodes eliminated
%  low_VOT_energy = zeros(np,nwin);
%  for iwin = 1:41
%      for ipat = 1:np
%          mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
%          mat = mat(current_ch_list,current_ch_list);
%          low_VOT_energy(ipat,iwin) = norm(mat,'fro');
%      end
%  end
%  
%  X1 = low_VOT_energy;
%  
%  stats = zeros(1,nwin);
%  tstat_val = zeros(1,nwin);
%  
%  ipat = np;
%  pat_list = 1:np;
%  pat_list(ipat)=[];
%  for iwin = 1:nwin
%      [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
%      stats(iwin) = p;
%      tstat_val(iwin) = t;
%  end
%  
%  %find the clusters, and max tstat cluster
%  tmask = stats<0.05;
%  cluster_starts = strfind((stats<0.05),[0 1])+1; %+1 as we dont want to include 0 in the string 011
%  if(tmask(1)==1)
%      cluster_starts = [1 cluster_starts];
%  end
%  
%  t_stats_clusters = zeros(1,length(cluster_starts));
%  for j = 1:length(cluster_starts)
%      within_cluster_incementer = cluster_starts(j);
%      while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
%          t_stats_clusters(j) = t_stats_clusters(j)+ t_statistics(within_cluster_incementer);
%          within_cluster_incementer = within_cluster_incementer+1;
%      end
%  end
%  
%  tstats_for_iterations(1) = sum(t_stats_clusters); % or max
% 
%  % now we start eliminating a node at a time.
%  
%  for ind_elim_nodes = 2:19
%      
%      tstats_for_sets_nodes = zeros(1,length(current_ch_list));
%      for inode = 1:length(current_ch_list)
%          
%          % cycle through the nodes to see which node to eliminate
%          temp_ch_list = current_ch_list; %all the nodes in consideration in this list
%          temp_ch_list(inode)=[]; % eliminate the current node.
%          low_VOT_energy = zeros(np,nwin);
%          for iwin = 1:41
%              for ipat = 1:np
%                  mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
%                  mat = mat(temp_ch_list,temp_ch_list);
%                  low_VOT_energy(ipat,iwin) = norm(mat,'fro');
%              end
%          end
%          
%          X1 = low_VOT_energy;
%          
%          stats = zeros(1,nwin);
%          tstat_val = zeros(1,nwin);
%          
%          ipat = np;
%          pat_list = 1:np;
%          pat_list(ipat)=[];
%          for iwin = 1:nwin
%              [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
%              stats(iwin) = p;
%              tstat_val(iwin) = t;
%          end
%          
%          %find the clusters, and max tstat cluster
%          tmask = stats<0.05;
%          cluster_starts = strfind((stats<0.05),[0 1])+1; %+1 as we dont want to include 0 in the string 011
%          if(tmask(1)==1)
%              cluster_starts = [1 cluster_starts];
%          end
%          
%          t_stats_clusters = zeros(1,length(cluster_starts));
%          for j = 1:length(cluster_starts)
%              within_cluster_incementer = cluster_starts(j);
%              while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
%                  t_stats_clusters(j) = t_stats_clusters(j)+ t_statistics(within_cluster_incementer);
%                  within_cluster_incementer = within_cluster_incementer+1;
%              end
%          end
%          
%          tstats_for_sets_nodes(inode) = sum(t_stats_clusters);
%      end
%      
%      [mm,ind] = max(tstats_for_sets_nodes);
%      
%      % eliminiate index ind form the current_ch_list
%      current_ch_list(ind)=[];
%      % save the max t-stat from this cluster
%      tstats_for_iterations(ind_elim_nodes) = mm;
%  end
%  
% tstats_for_iterations
% current_ch_list
% 
% 
% %% (2) new algorithm, delete a node at a time aphasia vs rest - lowVOT use max number of windows AND maxtstat
% 
% load('ChanNamesandLocs.mat');
% 
% k_n = 26;
% %X = cell2mat(chanNamesLocs(1,2:4));
% X = cell2mat(chanNamesLocs([1:25 28],2:4));
%     
% % find energy for each subplot, lowVOT, highVOT and midVOT
% current_ch_list = 1:26;
% tstats_for_iterations = zeros(1,length(chlist));
% 
% % full netowork stats: no nodes eliminated
%  low_VOT_energy = zeros(np,nwin);
%  for iwin = 1:41
%      for ipat = 1:np
%          mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
%          mat = mat(current_ch_list,current_ch_list);
%          low_VOT_energy(ipat,iwin) = norm(mat,'fro');
%      end
%  end
%  
%  X1 = low_VOT_energy;
%  
%  stats = zeros(1,nwin);
%  tstat_val = zeros(1,nwin);
%  
%  ipat = np;
%  pat_list = 1:np;
%  pat_list(ipat)=[];
%  for iwin = 1:nwin
%      [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
%      stats(iwin) = p;
%      tstat_val(iwin) = t;
%  end
%  
%  %find the clusters, and max tstat cluster
%  stats = fdr(stats);
%  tmask = stats<0.05;
%  cluster_starts = strfind((stats<0.05),[0 1])+1; %+1 as we dont want to include 0 in the string 011
%  if(tmask(1)==1)
%      cluster_starts = [1 cluster_starts];
%  end
%  
%  t_stats_clusters = zeros(1,length(cluster_starts));
%  for j = 1:length(cluster_starts)
%      within_cluster_incementer = cluster_starts(j);
%      while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
%          t_stats_clusters(j) = t_stats_clusters(j)+ t_statistics(within_cluster_incementer);
%          within_cluster_incementer = within_cluster_incementer+1;
%      end
%  end
%  
%  tstats_for_iterations(1) = max(t_stats_clusters).*sum(tmask); % or max
% 
%  % now we start eliminating a node at a time.
%  
%  for ind_elim_nodes = 2:22
%      
%      tstats_for_sets_nodes = zeros(1,length(current_ch_list));
%      for inode = 1:length(current_ch_list)
%          
%          % cycle through the nodes to see which node to eliminate
%          temp_ch_list = current_ch_list; %all the nodes in consideration in this list
%          temp_ch_list(inode)=[]; % eliminate the current node.
%          low_VOT_energy = zeros(np,nwin);
%          for iwin = 1:41
%              for ipat = 1:np
%                  mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
%                  mat = mat(temp_ch_list,temp_ch_list);
%                  low_VOT_energy(ipat,iwin) = norm(mat,'fro');
%              end
%          end
%          
%          X1 = low_VOT_energy;
%          
%          stats = zeros(1,nwin);
%          tstat_val = zeros(1,nwin);
%          
%          ipat = np;
%          pat_list = 1:np;
%          pat_list(ipat)=[];
%          for iwin = 1:nwin
%              [h,p, t, df] =  ttestch(mean(X1(pat_list,iwin)), std(X1(pat_list,iwin)), X1(ipat,iwin), np-1, 0.05);
%              stats(iwin) = p;
%              tstat_val(iwin) = t;
%          end
%          
%          %find the clusters, and max tstat cluster
%           stats = fdr(stats);
%          tmask = stats<0.05;
%          cluster_starts = strfind((stats<0.05),[0 1])+1; %+1 as we dont want to include 0 in the string 011
%          if(tmask(1)==1)
%              cluster_starts = [1 cluster_starts];
%          end
%          
%          t_stats_clusters = zeros(1,length(cluster_starts));
%          for j = 1:length(cluster_starts)
%              within_cluster_incementer = cluster_starts(j);
%              while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
%                  t_stats_clusters(j) = t_stats_clusters(j)+ t_statistics(within_cluster_incementer);
%                  within_cluster_incementer = within_cluster_incementer+1;
%              end
%          end
%          
%          tstats_for_sets_nodes(inode) = max(t_stats_clusters)*sum(tmask);
%      end
%      
%      [mm,ind] = max(tstats_for_sets_nodes);
%      
%      % eliminiate index ind form the current_ch_list
%      current_ch_list(ind)=[];
%      % save the max t-stat from this cluster
%      tstats_for_iterations(ind_elim_nodes) = mm;
%  end
%  
% tstats_for_iterations
% current_ch_list

%% specific brain regions lowVOT aphasic vs healthy

%X = cell2mat(chanNamesLocs(1,2:4));
%X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 4;

nearestNeighbors = cell(1,nbrainRegions);

  nearestNeighbors{1} = [1 2 3 4 5 6 7 9 10];
 nearestNeighbors{2} =     [ 11 15 16 19 20];
 nearestNeighbors{3} =     [8 13 17 18 12 14];
  nearestNeighbors{4} =     [21 22 23 24 25 26];
    

% find energy for each subplot

sig_windows_aphasia = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
    
    low_VOT_energy = zeros(np,nwin);
    for iwin = 1:41
        for ipat = 1:np
            mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
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
    str = sprintf('brain region %d',inode);
    title(str);
    box off;
    
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

    imagesc(ntim, 1:22,stats_low_gt_bs);
    sig_windows_aphasia(inode,:)=stats_low_gt_bs(end,:);
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
     xlabel('Time Windows');
     ylabel('subnetwork centered around node' );
     title('Aphisic vs healthy - using LowVOT measure');
     
%% specific brain regions highVOT aphasic vs healthy

%X = cell2mat(chanNamesLocs(1,2:4));
%X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 4;

nearestNeighbors = cell(1,nbrainRegions);

  nearestNeighbors{1} = [1 2 3 4 5 6 7 9 10];
 nearestNeighbors{2} =     [ 11 15 16 19 20];
 nearestNeighbors{3} =     [8 13 17 18 12 14];
  nearestNeighbors{4} =     [21 22 23 24 25 26];
    

% find energy for each subplot

sig_windows_aphasia = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
    
    low_VOT_energy = zeros(np,nwin);
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
    str = sprintf('brain region %d',inode);
    title(str);
    
    subplot(2,1,2);
    X1 = high_VOT_energy;

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

    imagesc(ntim, 1:22,stats_low_gt_bs);
    sig_windows_aphasia(inode,:)=stats_low_gt_bs(end,:);
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
     xlabel('Time Windows');
     ylabel('subnetwork centered around node' );
     title('Aphisic vs healthy - using LowVOT measure');
     
     

%% let's look at knn subgraphs, for significance between mid, low and high 

load('ChanNamesandLocs.mat');

%X = cell2mat(chanNamesLocs(1,2:4));
X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 4;

nearestNeighbors = cell(1,nbrainRegions);

  nearestNeighbors{1} = [1 2 3 4 5 6 7 9 10];
 nearestNeighbors{2} =     [ 11 15 16 19 20];
 nearestNeighbors{3} =     [8 13 17 18 12 14];
  nearestNeighbors{4} =     [21 22 23 24 25 26];
    
    
  
% find energy for each subplot, lowVOT, highVOT and midVOT

sig_windows_all = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
    
    high_VOT_energy_knn = zeros(np-1,nwin);
    low_VOT_energy_knn = zeros(np-1,nwin);
    mid_VOT_energy_knn = zeros(np-1,nwin);

    for iwin = 1:41
        for ipat = 1:np-1
            mat = squeeze(diff_highVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            high_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            low_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_midVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            mid_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
        end
    end
    
    pval_mid_gt_low_and_high = zeros(1,nwin);
    col = zeros(np-1,nwin);
    tstats = zeros(1,nwin);
    for i = 1:nwin
        anovaMatrix = [low_VOT_energy_knn(:,i) mid_VOT_energy_knn(:,i)  high_VOT_energy_knn(:,i)];
        bb = normalize_multiple_patients(anovaMatrix);
        col(:,i) = 2*bb(:,2) - (bb(:,1)+bb(:,3));
        [h,p,ci, stats] = ttest(col(:,i),0,'Tail','right');
        pval_mid_gt_low_and_high(i)=p;
        tstats(i) = stats.tstat;
    end
   % pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

    sig_windows_all(inode,:) = pval_mid_gt_low_and_high;
    
         %find the clusters, and max tstat cluster
         tmask = pval_mid_gt_low_and_high<0.05;
         cluster_starts = strfind((pval_mid_gt_low_and_high<0.05),[0 1])+1; %+1 as we dont want to include 0 in the string 011
         if(tmask(1)==1)
             cluster_starts = [1 cluster_starts];
         end
         
         t_stats_clusters = zeros(1,length(cluster_starts));
         for j = 1:length(cluster_starts)
             within_cluster_incementer = cluster_starts(j);
             while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
                 t_stats_clusters(j) = t_stats_clusters(j)+ tstats(within_cluster_incementer);
                 within_cluster_incementer = within_cluster_incementer+1;
             end
         end
         
         [mm,ind]=max(t_stats_clusters);
         max_cluster_start = cluster_starts(ind);
         
            nperm=500;
         [tstats_sorted,thresh] = func_cluster_analysis(nperm,nwin,np,low_VOT_energy_knn,mid_VOT_energy_knn,high_VOT_energy_knn,mm)

         cluster_end = strfind((pval_mid_gt_low_and_high<0.05),[1 0]); %+1 as we dont want to include 0 in the string 011
         if(tmask(end)==1)
             cluster_end = [cluster_end length(tmask)];
         end
         
        sig_cluster_starts = cluster_starts(t_stats_clusters>thresh);
        sig_cluster_end = cluster_end(t_stats_clusters>thresh);
        sig_pvalues = zeros(1,length(t_stats_clusters));
        for isig = 1:length(t_stats_clusters)            
            sig_pvalues(isig) = sum(tstats_sorted>t_stats_clusters(isig))./nperm;
        end
        sig_pvalues = sig_pvalues(t_stats_clusters>thresh);

    figure;
    plot(ntim,mean(col),'-*k','LineWidth',4);hold on;
    jbfill(ntim,mean(col)+std(col)*2/sqrt(21),mean(col)-std(col)*2/sqrt(21),'k','k',0,0.3);
    hold on;
    line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
    yl = ylim;

    
    for isig = 1:length(sig_cluster_starts)
        patch([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_starts(isig))], [yl(1) yl(1) yl(2) yl(2)],'y')
        alpha(0.3);
        str = sprintf('p = %0.2g',sig_pvalues(isig));
        text(mean([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig))]),yl(2)-0.01,str);
    end
end



%% mid > low and high? WHOLE brain


load('ChanNamesandLocs.mat');

%X = cell2mat(chanNamesLocs(1,2:4));
X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 1;

nearestNeighbors = cell(1,nbrainRegions);

  nearestNeighbors{1} = [1:26];

    
    
  
% find energy for each subplot, lowVOT, highVOT and midVOT

sig_windows_all = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
    
    high_VOT_energy_knn = zeros(np-1,nwin);
    low_VOT_energy_knn = zeros(np-1,nwin);
    mid_VOT_energy_knn = zeros(np-1,nwin);

    for iwin = 1:41
        for ipat = 1:np-1
            mat = squeeze(diff_highVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            high_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            low_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_midVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            mid_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
        end
    end
    
    pval_mid_gt_low_and_high = zeros(1,nwin);
    col = zeros(np-1,nwin);
    tstats = zeros(1,nwin);
    for i = 1:nwin
        anovaMatrix = [low_VOT_energy_knn(:,i) mid_VOT_energy_knn(:,i)  high_VOT_energy_knn(:,i)];
        bb = normalize_multiple_patients(anovaMatrix);
        col(:,i) = 2*bb(:,2) - (bb(:,1)+bb(:,3));
        [h,p,ci, stats] = ttest(col(:,i),0,'Tail','right');
        pval_mid_gt_low_and_high(i)=p;
        tstats(i) = stats.tstat;
    end
   % pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

    sig_windows_all(inode,:) = pval_mid_gt_low_and_high;
    
         %find the clusters, and max tstat cluster
         tmask = pval_mid_gt_low_and_high<0.05;
         cluster_starts = strfind((pval_mid_gt_low_and_high<0.05),[0 1])+1; %+1 as we dont want to include 0 in the string 011
         if(tmask(1)==1)
             cluster_starts = [1 cluster_starts];
         end
         
         t_stats_clusters = zeros(1,length(cluster_starts));
         for j = 1:length(cluster_starts)
             within_cluster_incementer = cluster_starts(j);
             while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
                 t_stats_clusters(j) = t_stats_clusters(j)+ tstats(within_cluster_incementer);
                 within_cluster_incementer = within_cluster_incementer+1;
             end
         end
         
         [mm,ind]=max(t_stats_clusters);
         max_cluster_start = cluster_starts(ind);
         
            nperm=500;
         [tstats_sorted,thresh] = func_cluster_analysis(nperm,nwin,np,low_VOT_energy_knn,mid_VOT_energy_knn,high_VOT_energy_knn,mm)

         cluster_end = strfind((pval_mid_gt_low_and_high<0.05),[1 0]); %+1 as we dont want to include 0 in the string 011
         if(tmask(end)==1)
             cluster_end = [cluster_end length(tmask)];
         end
         
        sig_cluster_starts = cluster_starts(t_stats_clusters>thresh);
        sig_cluster_end = cluster_end(t_stats_clusters>thresh);
        sig_pvalues = zeros(1,length(t_stats_clusters));
        for isig = 1:length(t_stats_clusters)            
            sig_pvalues(isig) = sum(tstats_sorted>t_stats_clusters(isig))./nperm;
        end
        sig_pvalues = sig_pvalues(t_stats_clusters>thresh);

    figure;
    plot(ntim,mean(col),'-*k','LineWidth',4);hold on;
    jbfill(ntim,mean(col)+std(col)*2/sqrt(21),mean(col)-std(col)*2/sqrt(21),'k','k',0,0.3);
    hold on;
    line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
    yl = ylim;

    
    for isig = 1:length(sig_cluster_starts)
        patch([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_starts(isig))], [yl(1) yl(1) yl(2) yl(2)],'y')
        alpha(0.3);
        str = sprintf('p = %0.3g',sig_pvalues(isig));
        text(mean([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig))]),yl(2)-0.01,str);
    end
end


%%  spectral clustering mid vs high and low

%func_cluster_connectivity_matrix(con_matrix, str,ch_names,k)
healthy_matricesLOW = squeeze(sum(diff_lowVOT(1:np-1,:,:,:)));
healthy_matricesHIGH = squeeze(sum(diff_highVOT(1:np-1,:,:,:)));
healthy_matricesMID = squeeze(sum(diff_midVOT(1:np-1,:,:,:)));
new_edge_matrices = healthy_matricesMID - 0.5*(healthy_matricesLOW+ healthy_matricesHIGH);
new_edge_matrices(new_edge_matrices<0)=0;


M.matrices = new_edge_matrices;

mat = (squeeze(M.matrices(:,:,1)));
n = numel(mat(:,1));
gamma = 0.9;
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

nbrainRegions = 3;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = find(ch==1);
nearestNeighbors{2} = find(ch==2);
nearestNeighbors{3} = find(ch==3);
% nearestNeighbors{4} = find(ch==4);
% nearestNeighbors{5} = find(ch==5);
% nearestNeighbors{6} = find(ch==6);
% nearestNeighbors{7} = find(ch==7);


sig_windows_aphasia = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
    
     high_VOT_energy_knn = zeros(np-1,nwin);
    low_VOT_energy_knn = zeros(np-1,nwin);
    mid_VOT_energy_knn = zeros(np-1,nwin);

    for iwin = 1:41
        for ipat = 1:np-1
            mat = squeeze(diff_highVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            high_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            low_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
            mat = squeeze(diff_midVOT(ipat,:,:,iwin));
            mat = mat(nearestNeighbors{inode},nearestNeighbors{inode});
            mid_VOT_energy_knn(ipat,iwin) = norm(mat,'fro');
            
        end
    end
        
    pval_mid_gt_low_and_high = zeros(1,nwin);
    col = zeros(np-1,nwin);
    tstats = zeros(1,nwin);
    for i = 1:nwin
        anovaMatrix = [low_VOT_energy_knn(:,i) mid_VOT_energy_knn(:,i)  high_VOT_energy_knn(:,i)];
        bb = normalize_multiple_patients(anovaMatrix);
        col(:,i) = 2*bb(:,2) - (bb(:,1)+bb(:,3));
        [h,p,ci, stats] = ttest(col(:,i),0,'Tail','right');
        pval_mid_gt_low_and_high(i)=p;
        tstats(i) = stats.tstat;
    end
   % pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

    sig_windows_all(inode,:) = pval_mid_gt_low_and_high;
    
         %find the clusters, and max tstat cluster
         tmask = pval_mid_gt_low_and_high<0.05;
         cluster_starts = strfind((pval_mid_gt_low_and_high<0.05),[0 1])+1; %+1 as we dont want to include 0 in the string 011
         if(tmask(1)==1)
             cluster_starts = [1 cluster_starts];
         end
         
         t_stats_clusters = zeros(1,length(cluster_starts));
         for j = 1:length(cluster_starts)
             within_cluster_incementer = cluster_starts(j);
             while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
                 t_stats_clusters(j) = t_stats_clusters(j)+ tstats(within_cluster_incementer);
                 within_cluster_incementer = within_cluster_incementer+1;
             end
         end
         
         [mm,ind]=max(t_stats_clusters);
         max_cluster_start = cluster_starts(ind);
         
         nperm=500;
         [tstats_sorted,thresh] = func_cluster_analysis(nperm,nwin,np,low_VOT_energy_knn,mid_VOT_energy_knn,high_VOT_energy_knn,mm)

         cluster_end = strfind((pval_mid_gt_low_and_high<0.05),[1 0]); %+1 as we dont want to include 0 in the string 011
         if(tmask(end)==1)
             cluster_end = [cluster_end length(tmask)];
         end
         
        sig_cluster_starts = cluster_starts(t_stats_clusters>thresh);
        sig_cluster_end = cluster_end(t_stats_clusters>thresh);
        sig_pvalues = zeros(1,length(t_stats_clusters));
        for isig = 1:length(t_stats_clusters)            
            sig_pvalues(isig) = sum(tstats_sorted>t_stats_clusters(isig))./nperm;
        end
        sig_pvalues = sig_pvalues(t_stats_clusters>thresh);

    figure;
    plot(ntim,mean(col),'-*k','LineWidth',4);hold on;
    jbfill(ntim,mean(col)+std(col)*2/sqrt(21),mean(col)-std(col)*2/sqrt(21),'k','k',0,0.3);
    hold on;
    line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
    yl = ylim;

    
    for isig = 1:length(sig_cluster_starts)
        patch([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_starts(isig))], [yl(1) yl(1) yl(2) yl(2)],'y')
        alpha(0.3);
        str = sprintf('p = %0.3g',sig_pvalues(isig));
        text(mean([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig))]),yl(2)-0.01,str);
    end
    
end

     %% plot this in the graphs

load('ChanNamesandLocs.mat');
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);
inode = 2;
chl = nearestNeighbors{inode};

empty_graph = zeros(length(chlist),length(chlist));
empty_graph1 = zeros(length(chl),length(chl));

% load the mat file
M.matrices = new_edge_matrices;
mat = (squeeze(M.matrices(:,:,1)));

empty_graph(chl,chl) = mat(chl,chl);

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

%% plot the empty graph for figure
chlist = [1:25 28];

load('ChanNamesandLocs.mat');
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);

empty_graph = zeros(length(chlist),length(chlist));


G = graph(empty_graph,chnames);
gf=figure;
p=plot(G,'XData',cell2mat(chanNamesLocs(chlist,2)),'YData',cell2mat(chanNamesLocs(chlist,3)),'NodeFontSize',14); hold on;
p.NodeColor='r';
hold on;
az=-90;
ele = 90;
view([az ele])
box off;
axis off;
set(gcf,'color','w');
annotation(gf,'ellipse',...
    [0.176,0.166666666666667,0.722214285714286,0.721428571428572]);

annotation(gf,'line',[0.473214285714286 0.516071428571428],...
    [0.880952380952381 0.933333333333334]);

annotation(gf,'line',[0.558928571428571 0.516071428571428],...
    [0.885714285714286 0.933333333333334]);



