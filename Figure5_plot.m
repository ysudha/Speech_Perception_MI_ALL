% spectral clusteromg for healthy midVOT- 0.5*(low + high)
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
% example of a partition with time
figure; imagesc('XData',[ntim],'YData',1:26,'CData',partition_multi);
ylim([1,26]);xlim([0,400]);
box off;
set(gcf,'color','w');
xlabel('Time (ms)','FontSize',16);
ylabel('Channel number','FontSize',16);
%colorbar('Ticks',[1,2,3]);
ax=gca;
ax.FontSize = 16;


ch = partition_multi(:,1);

%% let's repeat this partitioning multiple times, and save them

nperm = 500;
partition_indices = zeros(length(chlist),nperm);
for iperm = 1:nperm
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
    partition_indices(:,iperm) = partition_multi(:,1);    
end

%% consensus clustering
figure('Position',[100 100 1500 500]); 
imagesc(partition_indices);
box off;
set(gcf,'color','w');
xlabel('Repetitions','FontSize',16);
ylabel('Channel number','FontSize',16);
colorbar('Ticks',[1,2,3,4]);
ax=gca;
ax.FontSize = 16;

D= agreement(partition_indices);
figure; imagesc(D)
ciu = consensus_und(D,1,1000);
ch =ciu;
list_ch = 1:length(chlist);
[~,ind] = sort(ch);

figure; imagesc(D(list_ch(ind),list_ch(ind)))
box off;
set(gcf,'color','w');
xlabel('Channels','FontSize',16);
ylabel('Channels','FontSize',16);
axis off;
axis square;
colorbar;
ax=gca;
ax.FontSize = 16;
%% lets plot the energies for these individual partitions

ch=ciu;
nbrainRegions = 3;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = find(ch==1);
nearestNeighbors{2} = find(ch==2);
nearestNeighbors{3} = find(ch==3);

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
    legend({'MidVOT - 0.5(LowVOT + HighVOT)'});

    a=jbfill(ntim,mean(col)+std(col)*2/sqrt(21),mean(col)-std(col)*2/sqrt(21),'k','k',0,0.3);
    hold on;
    b=line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
    yl = ylim;
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    legend box off;
    set(gcf,'color','w');
    xlabel('Time (ms)','FontSize',16);
    ylabel('Contrast vector','FontSize',16);
    ax=gca;
    ax.FontSize = 16;
    box off;
    
    for isig = 1:length(sig_cluster_starts)
        a=patch([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_starts(isig))], [yl(1) yl(1) yl(2) yl(2)],'y')
        alpha(0.3);
        str = sprintf('p = %0.3g',sig_pvalues(isig));
        %text(mean([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig))]),yl(2)-0.01,str);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';

    end
    
    
    
        figure;
    plot(ntim,mean(col),'-*k','LineWidth',4);hold on;
    legend({'MidVOT - 0.5(LowVOT + HighVOT)'});

    a=jbfill(ntim,mean(col)+std(col)*2/sqrt(21),mean(col)-std(col)*2/sqrt(21),'k','k',0,0.3);
    hold on;
    b=line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
    yl = ylim;
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    legend box off;
    set(gcf,'color','w');
    xlabel('Time (ms)','FontSize',16);
    ylabel('Contrast vector','FontSize',16);
    ax=gca;
    ax.FontSize = 16;
    box off;
    
    for isig = 1:length(sig_cluster_starts)
        a=patch([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_end(isig)) ntim(sig_cluster_starts(isig))], [yl(1) yl(1) yl(2) yl(2)],'y')
        alpha(0.3);
        str = sprintf('p = %0.3g',sig_pvalues(isig));
        text(mean([ntim(sig_cluster_starts(isig)) ntim(sig_cluster_end(isig))]),yl(2)-0.01,str);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';

    end
    
    
end

     %% plot this in the graphs (cluster 1)

load('ChanNamesandLocs.mat');
chanNamesLocsShort = chanNamesLocs(chlist,:);
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);
inode = 1;
chl = nearestNeighbors{inode};

empty_graph = zeros(length(chlist),length(chlist));
graph1 = zeros(length(chl),length(chl));

% load the mat file
M.matrices = new_edge_matrices;

mat = (mean(M.matrices(:,:,:),3)); % mean edge metric

gf=figure;
empty_graph(chl,chl) = mat(chl,chl);

G = graph(empty_graph,chnames);
G1 = graph(graph1,chnames(chl));

LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
pe=plot(G,'XData',cell2mat(chanNamesLocsShort(:,2)),'YData',cell2mat(chanNamesLocsShort(:,3)),'NodeFontSize',14,'LineWidth',LWidths); hold on;
p=plot(G1,'XData',cell2mat(chanNamesLocsShort(chl,2)),'YData',cell2mat(chanNamesLocsShort(chl,3)),'NodeFontSize',14);


numlinks = degree(G);
numlinks(numlinks>0) = round(sum(mat(chl,chl)).*20);
highlight(pe,numlinks>0,'NodeColor','r')
pe.MarkerSize = [4+numlinks];

p.NodeColor='r';
az=-90;
ele = 90;
view([az ele])
box off;
set(gcf,'color','w');
annotation(gf,'ellipse',...
    [0.176,0.166666666666667,0.722214285714286,0.721428571428572]);

annotation(gf,'line',[0.473214285714286 0.516071428571428],...
    [0.880952380952381 0.933333333333334]);

annotation(gf,'line',[0.558928571428571 0.516071428571428],...
    [0.885714285714286 0.933333333333334]);
axis off;


     %% plot this in the graphs (cluster 2)

load('ChanNamesandLocs.mat');
chanNamesLocsShort = chanNamesLocs(chlist,:);
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);
inode = 2;
chl = nearestNeighbors{inode};

empty_graph = zeros(length(chlist),length(chlist));
graph1 = zeros(length(chl),length(chl));

% load the mat file
M.matrices = new_edge_matrices;

mat = (mean(M.matrices(:,:,:),3)); % mean edge metric

gf=figure;
empty_graph(chl,chl) = mat(chl,chl);

G = graph(empty_graph,chnames);
G1 = graph(graph1,chnames(chl));

LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
pe=plot(G,'XData',cell2mat(chanNamesLocsShort(:,2)),'YData',cell2mat(chanNamesLocsShort(:,3)),'NodeFontSize',14,'LineWidth',LWidths); hold on;
p=plot(G1,'XData',cell2mat(chanNamesLocsShort(chl,2)),'YData',cell2mat(chanNamesLocsShort(chl,3)),'NodeFontSize',14);


numlinks = degree(G);
numlinks(numlinks>0) = round(sum(mat(chl,chl)).*20);
highlight(pe,numlinks>0,'NodeColor','r')
pe.MarkerSize = [4+numlinks];

highlight(pe,numlinks>0,'NodeColor','r')
pe.MarkerSize = [4+numlinks];

p.NodeColor='r';
az=-90;
ele = 90;
view([az ele])
box off;
set(gcf,'color','w');
annotation(gf,'ellipse',...
    [0.176,0.166666666666667,0.722214285714286,0.721428571428572]);

annotation(gf,'line',[0.473214285714286 0.516071428571428],...
    [0.880952380952381 0.933333333333334]);

annotation(gf,'line',[0.558928571428571 0.516071428571428],...
    [0.885714285714286 0.933333333333334]);
axis off;


   %% plot this in the graphs (cluster 3)

load('ChanNamesandLocs.mat');
chanNamesLocsShort = chanNamesLocs(chlist,:);
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);
inode = 3;
chl = nearestNeighbors{inode};

empty_graph = zeros(length(chlist),length(chlist));
graph1 = zeros(length(chl),length(chl));

% load the mat file
M.matrices = new_edge_matrices;

mat = (mean(M.matrices(:,:,:),3)); % mean edge metric

gf=figure;
empty_graph(chl,chl) = mat(chl,chl);

G = graph(empty_graph,chnames);
G1 = graph(graph1,chnames(chl));

LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
pe=plot(G,'XData',cell2mat(chanNamesLocsShort(:,2)),'YData',cell2mat(chanNamesLocsShort(:,3)),'NodeFontSize',14,'LineWidth',LWidths); hold on;
p=plot(G1,'XData',cell2mat(chanNamesLocsShort(chl,2)),'YData',cell2mat(chanNamesLocsShort(chl,3)),'NodeFontSize',14);


numlinks = degree(G);
numlinks(numlinks>0) = round(sum(mat(chl,chl)).*20);
highlight(pe,numlinks>0,'NodeColor','r')
pe.MarkerSize = [4+numlinks];

p.NodeColor='r';
az=-90;
ele = 90;
view([az ele])
box off;
set(gcf,'color','w');
annotation(gf,'ellipse',...
    [0.176,0.166666666666667,0.722214285714286,0.721428571428572]);

annotation(gf,'line',[0.473214285714286 0.516071428571428],...
    [0.880952380952381 0.933333333333334]);

annotation(gf,'line',[0.558928571428571 0.516071428571428],...
    [0.885714285714286 0.933333333333334]);
axis off;
