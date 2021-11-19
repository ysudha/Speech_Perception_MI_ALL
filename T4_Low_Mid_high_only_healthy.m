% back to drawing board. How do the individual edges b,aseline normalized perform? Any
% differentiability?

clear;

 %% (1) load graphs, edge = lowVOT
chlist = [1:25 28];
miVOT123 = load('mi_matrix_VOT123_smallwin.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 13;
nch = length(chlist);
ntim = [[1:25:51] [101:25:326]];
ntim = ntim*2-202;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT
diff_lowVOT = zeros(np,nch,nch,nwin-3);

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
    for iwin = 4:nwin
       baseline = (mi123{1,1} + mi123{2,1})./2;
       mat2 = squeeze(mi123{iwin,1} - baseline); 
       mat = squeeze(mat2);
       diff_lowVOT(ipat,:,:,iwin-3) = mat(chlist,chlist);
    end

end

nwin = nwin-3;

diff_lowVOT(diff_lowVOT<0)=0; %< is the right one
diff_lowVOT = abs(diff_lowVOT);

low_VOT_energy = zeros(np,nwin);

for iwin = 1:nwin
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
    plot(ntim(4:end),low_VOT_energy(i,:),'b-*','LineWidth',4);
    hold on;
end
plot(ntim(4:end),low_VOT_energy(np,:),'r-*','LineWidth',4);
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
imagesc(ntim(4:end), 1:22,stats_low_gt_bs);
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
  %% (3) load graphs, edge = midVOT
chlist = [1:25 28];
miVOT456 = load('mi_matrix_VOT456_smallwin.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 13;
nch = length(chlist);
ntim = [[1:25:51] [101:25:326]];
ntim = ntim*2-202;
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
    for iwin = 4:nwin
       baseline = (mi456{1,1} + mi456{2,1})./2;
       mat2 = squeeze(mi456{iwin,1} - baseline); 
       mat = squeeze(mat2);
       diff_midVOT(ipat,:,:,iwin-3) = mat(chlist,chlist);
    end

end


nwin = nwin-3;
ntim = ntim(4:end);

diff_midVOT(diff_midVOT<0)=0; %< is the right one
diff_midVOT = abs(diff_midVOT);

mid_VOT_energy = zeros(np,nwin);

for iwin = 1:nwin
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
title('Significant values High-BS>0');
xlabel('Time');
ylabel('Subjects');
ax = gca;
ax.XGrid = 'on'
grid minor;
axis square;
colorbar;

%% let's look at individual indiffernces between low, mid and high for normals

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
for i = 1:nwin
    anovaMatrix = [low_VOT(:,i) mid_VOT(:,i)  high_VOT(:,i)];
   bb = normalize_multiple_patients(anovaMatrix);
   %bb = anovaMatrix;
    col = bb(:,2) - 0.5*(bb(:,1)+bb(:,3));
    contrasts(:,i) = col;
    [h,p] = ttest(col,0,'Tail','right');
    pval_mid_gt_low_and_high(i)=p;
end
figure;
plot(contrasts');hold on;
figure;
%plot(mean(contrasts),'k');
plot(ntim,mean(contrasts),'-*k','LineWidth',4);hold on;
jbfill(ntim,mean(contrasts)+std(contrasts)*2.086/sqrt(21),mean(contrasts)-std(contrasts)*2.086/sqrt(21),'k','k',0,0.3);
hold on;
line([0,ntim(end)],[0 0],'Color','blue','LineStyle','--');
xlabel('time'); ylabel('contrast vector');
figure;
bar(ntim,pval_mid_gt_low_and_high); hold on;
line([0,ntim(end)],[0.05 0.05],'Color','black','LineStyle','--');
ylabel('p-values');xlabel('starting point of time-windows')
%pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

%% comparing statistically at each time point, using ranova

low_VOT = healthy_low_VOT_energy;
high_VOT = healthy_high_VOT_energy;
mid_VOT = healthy_mid_VOT_energy;
pval_mid_gt_low_and_high = zeros(1,nwin);
for i = 1:nwin
    anovaMatrix = [low_VOT(:,i) mid_VOT(:,i)  high_VOT(:,i)];
   bb = normalize_multiple_patients(anovaMatrix);

   Ids =  cellstr(num2str([1:21]')); %id's of patients

    t = table(Ids,bb(:,1),bb(:,2),bb(:,3),'VariableNames',{'IDs','VOT1','VOT2','VOT3'});

Meas = dataset([1 2 3]','VarNames',{'Measurements'});

rm = fitrm(t,'VOT1-VOT3 ~ IDs-1','WithinDesign',Meas);

ranovatbl = ranova(rm);

rr = table2array(ranovatbl);

    pval_mid_gt_low_and_high(i)=fcdf(1./rr(1,4),rr(2,2),rr(1,2),'upper');
end

figure;
bar(pval_mid_gt_low_and_high)
%pval_mid_gt_low_and_high = fdr(pval_mid_gt_low_and_high);

%% let's look at knn subgraphs, for significance between mid, low and high 

load('ChanNamesandLocs.mat');

k_n = 15;
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


    for iwin = 1:nwin
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
     
%% comparing statistically only looking at means for all time

low_VOT = healthy_low_VOT_energy;
high_VOT = healthy_high_VOT_energy;
mid_VOT = healthy_mid_VOT_energy;

anovaMatrix = [mean(low_VOT(:,:),2) mean(mid_VOT(:,:),2)  mean(high_VOT(:,:),2)];
bb = normalize_multiple_patients(anovaMatrix);
[p,tbl, stats]=anova2(anovaMatrix)
c = multcompare(stats)
%col = bb(:,2) - 0.5*(bb(:,1)+bb(:,3));
%[h,p] = ttest(col);
p


%% let's look at knn subgraphs, for significance between mid, low and high 
% for mean of time points

load('ChanNamesandLocs.mat');

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
    

% find energy for each subplot, lowVOT, highVOT and midVOT

sig_windows_all = zeros(length(chlist),1);
    pval_mid_gt_low_and_high = zeros(1,length(chlist));

for inode = 1:length(chlist)
    
    high_VOT_energy_knn = zeros(np-1,nwin);
    low_VOT_energy_knn = zeros(np-1,nwin);
    mid_VOT_energy_knn = zeros(np-1,nwin);


    for iwin = 1:nwin
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
    sig_windows_all(inode) = p;

end


figure; %imagesc(sig_windows_aphasia)
[r, c] = size(sig_windows_all);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, sig_windows_all);          % Plot the image
%colormap(gray);                              % Use a gray colormap
%axis equal                                   % Make axes grid sizes equal
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');