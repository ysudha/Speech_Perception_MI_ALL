% plots healthy midVOT, highVOT, lowVOT.  good colors(fig 3)
% plots contrast vector, using cluster p-value (fig 3)


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
% figure;
% for i = 1:np-1
%     plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
%     hold on;
% end
% plot(ntim,low_VOT_energy(np,:),'r-*','LineWidth',4);
% ylabel('Frobenius norm');
% xlabel('Time');
% title('Energy of matrices: (LowVOT- BS) >0');


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

% figure;
% imagesc(ntim, 1:22,stats_low_gt_bs);
% title('Significant values Low-BS>0');
% xlabel('Time');
% ylabel('Subjects');
% ax = gca;
% ax.XGrid = 'on'
% grid minor;
% axis square;
% colorbar;

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
% figure;
% for i = 1:np-1
%     plot(ntim,high_VOT_energy(i,:),'b-*','LineWidth',4);
%     hold on;
% end
% plot(ntim,high_VOT_energy(np,:),'r-*','LineWidth',4);
% ylabel('Frobenius norm');
% xlabel('Time');
% title('Energy of matrices: (HighVOT- BS) >0');


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

% figure;
% imagesc(ntim, 1:22,stats_high_gt_bs);
% title('Significant values High-BS>0');
% xlabel('Time');
% ylabel('Subjects');
% ax = gca;
% ax.XGrid = 'on'
% grid minor;
% axis square;
% colorbar;
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
% figure;
% for i = 1:np-1
%     plot(ntim,mid_VOT_energy(i,:),'b-*','LineWidth',4);
%     hold on;
% end
% plot(ntim,mid_VOT_energy(np,:),'r-*','LineWidth',4);
% ylabel('Frobenius norm');
% xlabel('Time');
% title('Energy of matrices: (MidVOT- BS) >0');


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

% figure;
% imagesc(ntim, 1:22,stats_mid_gt_bs);
% title('Significant values Mid-BS>0');
% xlabel('Time');
% ylabel('Subjects');
% ax = gca;
% ax.XGrid = 'on'
% grid minor;
% axis square;
% colorbar;


%% let's look at individual indiffernces between low, mid and high for normals
% mean of healthy VOTS, and mean with confidence interval
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

% figure; hold on;
% plot(ntim,(aphasia_low_VOT_energy),'-.r','LineWidth',6);hold on;
% plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
% jbfill(ntim,mean(healthy_low_VOT_energy)+std(healthy_low_VOT_energy)*2.086/sqrt(21),mean(healthy_low_VOT_energy)-std(healthy_low_VOT_energy)*2.086/sqrt(21),'r','r',0,0.3);
% hold on;
% plot(ntim,(aphasia_high_VOT_energy),'-.b','LineWidth',6);hold on;
% plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
% jbfill(ntim,mean(healthy_high_VOT_energy)+std(healthy_high_VOT_energy)*2.086/sqrt(21),mean(healthy_high_VOT_energy)-std(healthy_high_VOT_energy)*2.086/sqrt(21),'b','b',0,0.3);
% hold on;
% plot(ntim,(aphasia_mid_VOT_energy),'-.k','LineWidth',6);hold on;
% plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
% jbfill(ntim,mean(healthy_mid_VOT_energy)+std(healthy_mid_VOT_energy)*2.086/sqrt(21),mean(healthy_mid_VOT_energy)-std(healthy_mid_VOT_energy)*2.086/sqrt(21),'k','k',0,0.3);
% hold on;
% legend('Low','','High','','Mid','');
% 
% figure;
% hold on;
% ntim = ([101:5:301]-101)*2;
% plot(ntim,(aphasia_low_VOT_energy),'-.r','LineWidth',6);hold on;
% plot(ntim,(aphasia_high_VOT_energy),'-.b','LineWidth',6);hold on;
% plot(ntim,(aphasia_mid_VOT_energy),'-.k','LineWidth',6);hold on;
% 
% plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
% plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
% plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
% legend('Low','High','Mid');
% title('Energy of VOT condition (mean of all healthy people at each time) vs aphasic (dotted line)');
% ylabel('Mean energy of MI matrices');
% xlabel('Time (ms)');

% %% paper figure for plotting: just healthy mean VOTs with error bars
% figure; hold on;
% plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
% errorbar(ntim,mean(healthy_low_VOT_energy),std(healthy_low_VOT_energy)*1/sqrt(21),'r');
% %jbfill(ntim,mean(healthy_low_VOT_energy)+std(healthy_low_VOT_energy)*2.086/sqrt(21),mean(healthy_low_VOT_energy)-std(healthy_low_VOT_energy)*2.086/sqrt(21),'r','r',0,0.3);
% hold on;
% plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
% errorbar(ntim,mean(healthy_high_VOT_energy),std(healthy_high_VOT_energy)*1/sqrt(21),'b');
% %jbfill(ntim,mean(healthy_high_VOT_energy)+std(healthy_high_VOT_energy)*2.086/sqrt(21),mean(healthy_high_VOT_energy)-std(healthy_high_VOT_energy)*2.086/sqrt(21),'b','b',0,0.3);
% hold on;
% plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
% errorbar(ntim,mean(healthy_mid_VOT_energy),std(healthy_mid_VOT_energy)*1/sqrt(21),'k');
% %jbfill(ntim,mean(healthy_mid_VOT_energy)+std(healthy_mid_VOT_energy)*2.086/sqrt(21),mean(healthy_mid_VOT_energy)-std(healthy_mid_VOT_energy)*2.086/sqrt(21),'k','k',0,0.3);
% hold on;
% legend('Low','','High','','Mid','');


%% paper figure for plotting: just healthy mean VOTs with error bars
figure; hold on;
ntim = ([101:5:301]-101)*2;
ntim = ntim+50; % mean of window
plot(ntim,mean(healthy_low_VOT_energy),'-d','LineWidth',4,'MarkerSize',4,'Color',[27,158,119]/255);hold on;
plot(ntim,mean(healthy_high_VOT_energy),'-s','LineWidth',4,'MarkerSize',4,'Color',[117,112,179]/255); hold on;
hold on;
plot(ntim,mean(healthy_mid_VOT_energy),'-o','LineWidth',4,'MarkerSize',4,'Color',[217,95,2]/255);hold on;
hold on;
legend({'LowVOT','HighVOT','MidVOT'});

hold on;
a=jbfill(ntim,mean(healthy_high_VOT_energy)+std(healthy_high_VOT_energy)*1/sqrt(21),mean(healthy_high_VOT_energy)-std(healthy_high_VOT_energy)*1/sqrt(21),[117,112,179]/255,[117,112,179]/255,0,0.1);
hold on;
b=jbfill(ntim,mean(healthy_low_VOT_energy)+std(healthy_low_VOT_energy)*1/sqrt(21),mean(healthy_low_VOT_energy)-std(healthy_low_VOT_energy)*1/sqrt(21),[27,158,119]/255,[27,158,119]/255,0,0.1);
hold on;
c=jbfill(ntim,mean(healthy_mid_VOT_energy)+std(healthy_mid_VOT_energy)*1/sqrt(21),mean(healthy_mid_VOT_energy)-std(healthy_mid_VOT_energy)*1/sqrt(21),[217,95,2]/255,[217,95,2]/255,0,0.1);
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
b.Annotation.LegendInformation.IconDisplayStyle = 'off';
c.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlim([25,475]);
%ylim([0.24 0.48]);
box off;
set(gcf,'color','w');
xlabel('Time (ms)','FontSize',16);
ylabel('Average Frobenius norm for healthy subjects','FontSize',16);
ax=gca;
ax.FontSize = 16;
legend box off

%% mid > low and high? WHOLE brain


load('ChanNamesandLocs.mat');

%X = cell2mat(chanNamesLocs(1,2:4));
X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 1;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = [1:26];

% find energy for each subplot, lowVOT, highVOT and midVOT

sig_windows_all = zeros(nbrainRegions,nwin);
inode = 1;

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

%% plotting the figure of contrast vector, with significance testing
ntim = ([101:5:301]-101)*2;
ntim = ntim+50; % mean of window
figure;
plot(ntim,mean(col),'-*k','LineWidth',4);hold on;
xlim([25,475]);
legend({'MidVOT - 0.5(LowVOT + HighVOT)'});
a=jbfill(ntim,mean(col)+std(col)*1/sqrt(21),mean(col)-std(col)*1/sqrt(21),'k','k',0,0.2);
hold on;
b=line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
b.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend box off
yl = ylim;

box off;
set(gcf,'color','w');
xlabel('Time (ms)','FontSize',16);
ylabel('Contrast vector','FontSize',16);
ax=gca;
ax.FontSize = 16;

for isig = 1:length(sig_cluster_starts)
    a=patch([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_starts(isig))], [yl(1) yl(1) yl(2) yl(2)],'y')
    alpha(0.3);
    str = sprintf('**',sig_pvalues(isig));
    text(mean([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig))])-6,yl(2)+0.01,str,'FontSize',20);
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';

end
