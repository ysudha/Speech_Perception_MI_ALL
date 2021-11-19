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

%% paper figure for plotting: just healthy mean VOTs with error bars
figure; hold on;
plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
errorbar(ntim,mean(healthy_low_VOT_energy),std(healthy_low_VOT_energy)*1/sqrt(21),'r');
%jbfill(ntim,mean(healthy_low_VOT_energy)+std(healthy_low_VOT_energy)*2.086/sqrt(21),mean(healthy_low_VOT_energy)-std(healthy_low_VOT_energy)*2.086/sqrt(21),'r','r',0,0.3);
hold on;
plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
errorbar(ntim,mean(healthy_high_VOT_energy),std(healthy_high_VOT_energy)*1/sqrt(21),'b');
%jbfill(ntim,mean(healthy_high_VOT_energy)+std(healthy_high_VOT_energy)*2.086/sqrt(21),mean(healthy_high_VOT_energy)-std(healthy_high_VOT_energy)*2.086/sqrt(21),'b','b',0,0.3);
hold on;
plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
errorbar(ntim,mean(healthy_mid_VOT_energy),std(healthy_mid_VOT_energy)*1/sqrt(21),'k');
%jbfill(ntim,mean(healthy_mid_VOT_energy)+std(healthy_mid_VOT_energy)*2.086/sqrt(21),mean(healthy_mid_VOT_energy)-std(healthy_mid_VOT_energy)*2.086/sqrt(21),'k','k',0,0.3);
hold on;
legend('Low','','High','','Mid','');


%% paper figure for plotting: just healthy mean VOTs with error bars
figure; hold on;
plot(ntim,mean(healthy_low_VOT_energy),'-*r','LineWidth',4);hold on;
jbfill(ntim,mean(healthy_low_VOT_energy)+std(healthy_low_VOT_energy)*1/sqrt(21),mean(healthy_low_VOT_energy)-std(healthy_low_VOT_energy)*1/sqrt(21),'r','r',0,0.1);
hold on;
plot(ntim,mean(healthy_high_VOT_energy),'-*b','LineWidth',4); hold on;
jbfill(ntim,mean(healthy_high_VOT_energy)+std(healthy_high_VOT_energy)*1/sqrt(21),mean(healthy_high_VOT_energy)-std(healthy_high_VOT_energy)*1/sqrt(21),'b','b',0,0.1);
hold on;
plot(ntim,mean(healthy_mid_VOT_energy),'-*k','LineWidth',4);hold on;
jbfill(ntim,mean(healthy_mid_VOT_energy)+std(healthy_mid_VOT_energy)*1/sqrt(21),mean(healthy_mid_VOT_energy)-std(healthy_mid_VOT_energy)*1/sqrt(21),'k','k',0,0.1);
hold on;
legend('Low','','High','','Mid','');
ylim([0.25 0.47]);
box off;
set(gcf,'color','w');
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

nperm = 500;
max_tstats = zeros(1,nperm);
for iperm = 1:nperm
    pval_mid_gt_low_and_high = zeros(1,nwin);
    contrasts = zeros(np-1,nwin);
    t_statistics = zeros(1,nwin);
    for i = 1:nwin
        anovaMatrix = [low_VOT(:,i) mid_VOT(:,i)  high_VOT(:,i)];
        bb = normalize_multiple_patients(anovaMatrix);
        for j = 1:np-1
            bb(j,:) = bb(j,randperm(3));
        end
        %bb = anovaMatrix;
        col = bb(:,2) - 0.5*(bb(:,1)+bb(:,3));
        contrasts(:,i) = col;
        [h,p,ci,stats] = ttest(col,0,'Tail','right');
        pval_mid_gt_low_and_high(i)=p;
        t_statistics(i) = stats.tstat;
    end
%     figure;
%     plot(contrasts');hold on;
%     figure;
%     %plot(mean(contrasts),'k');
%     plot(ntim,mean(contrasts),'-*k','LineWidth',4);hold on;
%     jbfill(ntim,mean(contrasts)+std(contrasts)*2.086/sqrt(21),mean(contrasts)-std(contrasts)*2.086/sqrt(21),'k','k',0,0.3);
%     hold on;
%     line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
%     xlabel('time'); ylabel('contrast vector');
%     figure;
%     bar(ntim,pval_mid_gt_low_and_high); hold on;
%     line([0,ntim(end)],[0.05 0.05],'Color','black','LineStyle','--');
%     ylabel('p-values');xlabel('starting point of time-windows')
%     
%     figure;
%     bar(ntim,t_statistics); hold on;
%     line([0,ntim(end)],[1.76 1.76],'Color','black','LineStyle','--');
%     ylabel('t-statistic');xlabel('starting point of time-windows')
tmask = t_statistics>1.725;
cluster_starts = strfind((t_statistics>1.725),[0 1])+1; %+1 as we dont want to include 0 in the string 011
t_stats_clusters = zeros(1,length(cluster_starts));
for j = 1:length(cluster_starts)
    
    within_cluster_incementer = cluster_starts(j);
    while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
        t_stats_clusters(j) = t_stats_clusters(j)+ t_statistics(within_cluster_incementer);
        within_cluster_incementer = within_cluster_incementer+1;
    end
end
  
if(length(t_stats_clusters)>0)
    max_tstats(iperm) = max(t_stats_clusters);
else
    max_tstats(iperm) = 0;
end

end

figure;
(hist(max_tstats));
tstats_sorted = sort(max_tstats);
thresh = tstats_sorted(round(nperm*0.95))

pval = sum(max_tstats>max_t_cluster)./nperm

%% let's look at knn subgraphs, for significance between mid, low and high 

load('ChanNamesandLocs.mat');

k_n = 8;
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
    

% find energy for each subplot, lowVOT, highVOT and midVOT

sig_windows_all = zeros(length(chlist),nwin);
for inode = 1:length(chlist)
    
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
[r, c] = size(sig_windows_all);    
%c=2;(:,9:10)
% Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_all);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
     
figure;
aa = sig_windows_all<1;
[r, c] = size(sig_windows_all);    
%c=2;(:,9:10)
% Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, aa);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
     

%% let's look at knn subgraphs, for significance between mid, low and high 
% for mean of time points

load('ChanNamesandLocs.mat');

k_n = 8;
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
    

% find energy for each subplot, lowVOT, highVOT and midVOT

sig_windows_all = zeros(1,length(chlist));
    pval_mid_gt_low_and_high = zeros(1,length(chlist));

for inode = 1:length(chlist)
    
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
    
  anovaMatrix = [mean(low_VOT_energy_knn(:,:),2) mean(mid_VOT_energy_knn(:,:),2)  mean(high_VOT_energy_knn(:,:),2)];
bb = normalize_multiple_patients(anovaMatrix);
%[p,tbl, stats]=anova1(anovaMatrix)
col = bb(:,2) - 0.5*(bb(:,1)+bb(:,3));
[h,p] = ttest(col,0,'Tail','right');
    sig_windows_all(1,inode) = p;

end


figure; %imagesc(sig_windows_aphasia)
[r, c] = size(sig_windows_all);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_all);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
     
     
    
%% specific brain regions LowVOT aphasia vs healthy
load('ChanNamesandLocs.mat');

k_n = 26;
%X = cell2mat(chanNamesLocs(1,2:4));
X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 16;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = [3 4 8 9 13];
nearestNeighbors{2} = [12 14 17 18 23  ];
nearestNeighbors{3} = [10 11 15 16 20];
nearestNeighbors{4} = [21 22 26];
nearestNeighbors{5} = [1 2 5 6 7];
nearestNeighbors{6} = [19 24 25 26]; 
nearestNeighbors{7} = [8 13 18 17]; % next knn add 9, next one add 22
nearestNeighbors{8} = [8 9 13 18 17 22]; 
 nearestNeighbors{9} = [2     3     5     9    11    12    17    18    23];
 nearestNeighbors{10} =     [  6    20    21    22    24    25];
 nearestNeighbors{11} =     [13    14    16    19];
  nearestNeighbors{12} =     [     1     4     7     8    10    15    26];
 
    nearestNeighbors{13} = [1 2 3 4 5 6 7 9 10];
 nearestNeighbors{14} =     [ 11 15 16 19 20];
 nearestNeighbors{15} =     [8 13 17 18 12 14];
  nearestNeighbors{16} =     [21 22 23 24 25 26];
  % find energy for each subplot, lowVOT, highVOT and midVOT


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
 
end

figure; %imagesc(sig_windows_aphasia)
[r, c] = size(sig_windows_aphasia);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_aphasia);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

     
     %% specific brain regions highVOT aphasia vs healthy
load('ChanNamesandLocs.mat');

%X = cell2mat(chanNamesLocs(1,2:4));
X = cell2mat(chanNamesLocs([1:25 28],2:4));
nbrainRegions = 16;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = [3 4 8 9 13];
nearestNeighbors{2} = [12 14 17 18 23  ];
nearestNeighbors{3} = [10 11 15 16 20];
nearestNeighbors{4} = [21 22 26];
nearestNeighbors{5} = [1 2 5 6 7];
nearestNeighbors{6} = [19 24 25 26]; 
nearestNeighbors{7} = [8 13 18 17]; % next knn add 9, next one add 22
nearestNeighbors{8} = [8 9 13 18 17 22]; 
 nearestNeighbors{9} = [2     3     5     9    11    12    17    18    23];
 nearestNeighbors{10} =     [  6    20    21    22    24    25];
 nearestNeighbors{11} =     [13    14    16    19];
  nearestNeighbors{12} =     [     1     4     7     8    10    15    26];
 
    nearestNeighbors{13} = [1 2 3 4 5 6 7 9 10];
 nearestNeighbors{14} =     [ 11 15 16 19 20];
 nearestNeighbors{15} =     [8 13 17 18 12 14];
  nearestNeighbors{16} =     [21 22 23 24 25 26];
  % find energy for each subplot, lowVOT, highVOT and midVOT


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

    sig_windows_aphasia(inode,:)=stats_low_gt_bs(end,:);
 
end

figure; %imagesc(sig_windows_aphasia)
[r, c] = size(sig_windows_aphasia);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_aphasia);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
