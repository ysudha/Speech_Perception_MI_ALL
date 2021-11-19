% Low - BS >0 and High - BS>0 edges are good. let's see what robust pca can
% give us. (1) robust pca sparse matrices are similar to the regular
% matrices. (2) node degrees, (3) which edges have greater average strength
% than other patients?



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

% what's the energy of the individual matrices look like? Something unique
% about baseline?
energy_ofLowVOT = zeros(np,43);

for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    for iwin = 1:43
       mat = squeeze(mi123{iwin,1});
       mat = mat(chlist,chlist);
       energy_ofLowVOT(ipat,iwin) = norm(mat,'fro');
    end
end

figure;
for i = 1:np-1
    plot([-200 -100 ntim],energy_ofLowVOT(i,:),'b-*','LineWidth',4);
    hold on;
end
plot([-200 -100 ntim],energy_ofLowVOT(np,:),'r-*','LineWidth',4);
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: LowVOT ');

% 
% %Why? Leys look at the normalized matrices themselves
% for iwin = 10
%     figure;
%     for ipat = 1:np    
%         mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
%         subplot(5,5,ipat);
%         imagesc(mat);
%         colormap jet;
%         caxis([0 0.2]);
%         axis square;
%     end
% end
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


%% plot energy of sparse matrices
energy_S_low_gt_BS = zeros(np,nwin);
S_mat = cell(1,np);
Lmat = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_lowVOT(ipat,:,:,1:nwin));
    X(X<0)=0;
    [Lx,Sx] = RobustPCA_time(X);
    S_mat{1,ipat} = Sx;
    Lmat{1,ipat} = Lx;
    for iwin = 1:nwin
       energy_S_low_gt_BS(ipat,iwin) = norm(squeeze(Sx(:,:,iwin)),'fro');
    end
    toc
end


X1 = energy_S_low_gt_BS(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-BS >0');

%Why? Leys look at the normalized matrices themselves
for iwin = 1:nwin
    figure;
    for ipat = 1:np    
        mat = (S_mat{1,ipat});
        mat = squeeze(mat(:,:,iwin));
        subplot(5,5,ipat);
        imagesc(mat);
        colormap jet;colorbar;
        caxis([0 0.15]);
        axis square;
    end
end

figure;
for ipat = 1:np
    mat = (Lmat{1,ipat});
    subplot(5,5,ipat);
    imagesc(mat);
    colormap jet;
    caxis([0 0.05]);
    axis square;colorbar;
end

L_energy = zeros(1,np);
figure;hold on;
for ipat = 1:np
    mat = (Lmat{1,ipat});
    bar(ipat, norm(mat,'fro'));
    L_energy(ipat) = norm(mat,'fro');
end


pval = zeros(1,np);
for ipat = 1:np
[h,p, t, df] =  ttestch(mean(L_energy), std(L_energy), L_energy(ipat), np-1, 0.05);
fprintf('Pat %d p<0.05 = %d\n',ipat, p);
pval(ipat) = p;
end
[p_fdr, p_masked] = fdr( pval, 0.05);

%% statistic for sparse

X1 = energy_S_low_gt_BS;

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


%% let's isolate by brain region

% regions to be grouped
% [1,2] [ 3,4], [5,9,10,14], [6,7],
% [8,13,18],[11,15,19],[12],[16],[17,21],[20,25], [22,23,24,26,27,28,29,30]
grp = {[1,2], [ 3,4], [5,9,10,14], [6,7],[8,13,18],[11,15,19],[12],[16],[17,21],[20,25],[22,23,24,28]};%, [22,23,24,26,27,28,29,30]};

%inf_front = [1,3,4,8,12,13];
%chlist = [1 2 4 6 7 9 11 12 13 14 16 17 18 26];
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

% what's the energy of the individual matrices look like? Something unique
% about baseline?
energy_ofLowVOT = zeros(np,43);

for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    for iwin = 1:43
       mat = squeeze(mi123{iwin,1});
       mat = mat(chlist,chlist);
       energy_ofLowVOT(ipat,iwin) = norm(mat,'fro');
    end
end

figure;
for i = 1:np-1
    plot([-200 -100 ntim],energy_ofLowVOT(i,:),'b-*','LineWidth',4);
    hold on;
end
plot([-200 -100 ntim],energy_ofLowVOT(np,:),'r-*','LineWidth',4);
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: LowVOT ');

% plot degrees of all patients
% for ipat = 1:22
% figure; hold on;
% for i = 1:41;
% plot(sum(squeeze(diff_lowVOT(ipat,:,:,i))))
% end
% end

% need a bar plot, avegerage degree of each node, for all patients
deg_pat = zeros(ipat,length(chlist));

for ipat = 1:22
    deg = zeros(1,length(chlist));
    for iwin = 1:nwin
        deg = deg + (sum(squeeze(diff_lowVOT(ipat,:,:,iwin))));
    end
    deg = deg./nwin;
    deg_pat(ipat,:) = (deg);
 
end

figure;
imagesc(deg_pat)

%% statistics for node degree

X1 = deg_pat;
nnodes = length(chlist);
stats = zeros(1,nnodes);
stats_l_bs = zeros(1,np);
stats_low_gt_bs = zeros(np,nnodes);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for in = 1:nnodes
        [h,p, t, df] =  ttestch(mean(X1(pat_list,in)), std(X1(pat_list,in)), X1(ipat,in), np-1, 0.05);
        stats(in) = p;
        
    end
    stats = fdr(stats);
    fprintf('low_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_l_bs(ipat) = sum(stats<0.05);
    stats_low_gt_bs(ipat,:) = stats<0.05;
end

figure;
imagesc(1:nnodes, 1:22,stats_low_gt_bs);
title('Significant values Low-BS>0');
xlabel('Time');
ylabel('Subjects');
ax = gca;
ax.XGrid = 'on'
grid minor;
axis square;
colorbar;

%% individual edges


% regions to be grouped
% [1,2] [ 3,4], [5,9,10,14], [6,7],
% [8,13,18],[11,15,19],[12],[16],[17,21],[20,25], [22,23,24,26,27,28,29,30]
grp = {[1,2], [ 3,4], [5,9,10,14], [6,7],[8,13,18],[11,15,19],[12],[16],[17,21],[20,25],[22,23,24,28]};%, [22,23,24,26,27,28,29,30]};

%inf_front = [1,3,4,8,12,13];
%chlist = [1 2 4 6 7 9 11 12 13 14 16 17 18 26];
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

% need a bar plot, avegerage strength of each edge, for all patients
mm = triu(squeeze(diff_lowVOT(1,:,:,1)),1);
ll =  length(mm(:));
edge_pat = zeros(ipat,ll);
mat = zeros(size(mm));
for ipat = 1:22
    deg = zeros(1,ll);
    for iwin = 1:nwin
        temp = triu(squeeze(diff_lowVOT(ipat,:,:,iwin)),1);
        mat = mat+temp;
    end
    mat = mat./nwin;
    edge_pat(ipat,:) = mat(:);
 
end

vec = find(sum(edge_pat)>0);
figure;
imagesc(edge_pat(:,vec))


%% statistics for individual edges

X1 = edge_pat(:,vec);
nedges = length(vec);
stats = zeros(1,nedges);
stats_l_bs = zeros(1,np);
stats_low_gt_bs = zeros(np,nedges);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    test = sum(edge_pat);
    for in = 1:nedges
        if(test(in)>0)
            %[h,p, t, df] =  ttestch(mean(X1(pat_list,in)), std(X1(pat_list,in)), X1(ipat,in), np-1, 0.05);
            [h,p, t, df] =  ttestch(mean(X1(pat_list,in)), std(X1(pat_list,in)), X1(ipat,in), np-1, 0.05);

            stats(in) = p;
        end
        
    end
    stats = fdr(stats);
    fprintf('low_less_BS %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));
    stats_l_bs(ipat) = sum(stats<0.05);
    stats_low_gt_bs(ipat,:) = stats<0.05;
end

figure;
imagesc(1:nedges, 1:22,stats_low_gt_bs);
title('Significant values Low-BS>0');
xlabel('Time');
ylabel('Subjects');
ax = gca;
ax.XGrid = 'on'
grid minor;
axis square;
colorbar;

