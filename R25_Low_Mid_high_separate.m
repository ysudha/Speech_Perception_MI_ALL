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


%Why? Leys look at the normalized matrices themselves
for iwin = 10
    figure;
    for ipat = 1:np    
        mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
        subplot(5,5,ipat);
        imagesc(mat);
        colormap jet;
        caxis([0 0.2]);
        axis square;
    end
end


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

%% edges in the actual matrix
% find the edges with greatest strength
mask_all = zeros(nch,nch,nwin);
for iwin = 1:41
    %figure; hold on;
    mask = ones(nch,nch);
    for ipat = 1:np-1
        S = squeeze(diff_lowVOT(ipat,:,:,iwin));
        Sap = squeeze(diff_lowVOT(22,:,:,iwin));
        mask_subj_specific = Sap>S;
        mask = mask.*mask_subj_specific;
    end
    %imagesc(mask)
    mask_all(:,:,iwin) = mask;
end

% relevant time windows which were significant
sig_win = [1:15 19:41];
sig_win = [1:10];

mask = squeeze(mask_all(:,:,sig_win(1)));
for i = 2:length(sig_win)
    mask = mask.*squeeze(mask_all(:,:,sig_win(i)));
end
figure;
imagesc(mask)

% plot on graph

load('ChanNamesandLocs.mat');
chnames = chanNamesLocs(:,1);
chnames = chnames(chlist);

G = graph(mask,chnames);
plot(G,'XData',cell2mat(chanNamesLocs(chlist,2)),'YData',cell2mat(chanNamesLocs(chlist,3)),'ZData',cell2mat(chanNamesLocs(chlist,4)))

%% new thing to try: how do the connections in the mask region different for 
% aphasia guy versus the rest
for iwin = 10
    figure; hold on;
    mask = ones(nch,nch);
    for ipat = 1:np
        S = squeeze(diff_lowVOT(ipat,:,:,iwin));
        plot(sum(S));
    end
    %imagesc(mask)
end
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

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    for iwin = 3:nwin
       mat2 = squeeze(mi789{iwin,1} - mi789{1,1}); 
       mat = squeeze(mat2);
       diff_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

diff_highVOT(diff_highVOT<0)=0; %< is correct
diff_highVOT = abs(diff_highVOT);

high_VOT_energy = zeros(np,nwin);

for iwin = 1:41
    %figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_highVOT(ipat,:,:,iwin));
        %bar(ipat,norm(mat,'fro'));
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
%title('Energy of matrices: (HighVOT- BS) <0');

% what's the energy of the individual matrices look like? Something unique
% about baseline?
energy_ofHighVOT = zeros(np,43);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    for iwin = 1:43
       mat = squeeze(mi789{iwin,1});
       mat = mat(chlist,chlist);
       energy_ofHighVOT(ipat,iwin) = norm(mat,'fro');
    end
end

figure;
for i = 1:np-1
    plot([-200 -100 ntim],energy_ofHighVOT(i,:),'b-*','LineWidth',4);
    hold on;
end
plot([-200 -100 ntim],energy_ofHighVOT(np,:),'r-*','LineWidth',4);
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: HighVOT ');

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
miVOT456 = load('mi_matrix_VOT456_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
ntim = [101:5:301]*2-202;

% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT
diff_midVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1}; %contains 43 time windows
    for iwin = 3:nwin
       mat2 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat = squeeze(mat2);
       diff_midVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

diff_midVOT(diff_midVOT<0)=0;
diff_midVOT = abs(diff_midVOT);

mid_VOT_energy = zeros(np,nwin);

for iwin = 1:41
    %figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_midVOT(ipat,:,:,iwin));
        %bar(ipat,norm(mat,'fro'));
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
%title('Energy of matrices: (MidVOT- BS) <0');

% what's the energy of the individual matrices look like? Something unique
% about baseline?
energy_ofMidVOT = zeros(np,43);

for ipat = 1:np
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1}; %contains 43 time windows
    for iwin = 1:43
       mat = squeeze(mi456{iwin,1});
       mat = mat(chlist,chlist);
       energy_ofMidVOT(ipat,iwin) = norm(mat,'fro');
    end
end

figure;
for i = 1:np-1
    plot([-200 -100 ntim],energy_ofMidVOT(i,:),'b-*','LineWidth',4);
    hold on;
end
plot([-200 -100 ntim],energy_ofMidVOT(np,:),'r-*','LineWidth',4);
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: MidVOT ');
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










