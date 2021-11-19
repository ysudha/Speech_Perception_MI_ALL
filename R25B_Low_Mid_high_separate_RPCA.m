% for asilomar paper
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

energy_Stime = zeros(np,nwin);
energy_Xtime = zeros(np,nwin);
S_pat_l = zeros(np,nch,nch,nwin);
L_pat_l = zeros(np,nch,nch);
parfor ipat = 1:np
    tic
    X = squeeze(diff_lowVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    S_pat_l(ipat,:,:,:)=Sx;
    L_pat_l(ipat,:,:)=Lx;
    toc
end



low_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(S_pat_l(ipat,:,:,iwin));
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
%title('Energy of matrices: LowVOT condition');
%title('Energy of matrices: (LowVOT- BS) <0');
set(gcf,'color','w');
box off;
ax=gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

figure;
for ipat = 1:np
     subplot(5,5,ipat);
     imagesc(squeeze(L_pat_l(ipat,:,:)));
     caxis([0 0.03]);
     axis square;
end


figure('position',[440,378,900,420]);
%figure;
L_energy = zeros(1,np);
for ipat = 1:np
    temp = squeeze(L_pat_l(ipat,:,:));
    temp(temp<0)=0;
    L_energy(ipat) = norm(temp,'fro');
end
h = bar(1:21,L_energy(1:21),'b');
hold on;
h1 = bar(22,L_energy(22),'r');
set(gcf,'color','w');
box off;
ax=gca;
ax.XTick = 1:22;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;



% statistics
[h,p, t, df] =  ttestch(mean(L_energy), std(L_energy), L_energy(end), np-1, 0.05)
p

%% (2) load graphs, edge = highVOT
chlist = [1:25 28];
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
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

diff_highVOT(diff_highVOT<0)=0; %< is the right one
diff_highVOT = abs(diff_highVOT);

energy_Stime = zeros(np,nwin);
energy_Xtime = zeros(np,nwin);
S_pat_h = zeros(np,nch,nch,nwin);
L_pat_h = zeros(np,nch,nch);
parfor ipat = 1:np
    tic
    X = squeeze(diff_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    S_pat_h(ipat,:,:,:)=Sx;
    L_pat_h(ipat,:,:)=Lx;
    toc
end



high_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(S_pat_h(ipat,:,:,iwin));
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
%title('Energy of matrices: HighVOT condition');
set(gcf,'color','w');
box off;
ax=gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;


figure;
for ipat = 1:np
     subplot(5,5,ipat);
     imagesc(squeeze(L_pat_h(ipat,:,:)));
     caxis([0 0.03]);

end

figure('position',[440,378,900,420]);
%figure;
L_energy = zeros(1,np);
for ipat = 1:np
    temp = squeeze(L_pat_h(ipat,:,:));
    temp(temp<0)=0;
    L_energy(ipat) = norm(temp,'fro');
end
h = bar(1:21,L_energy(1:21),'b');
hold on;
h1 = bar(22,L_energy(22),'r');
set(gcf,'color','w');
box off;
ax=gca;
ax.XTick = 1:22;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

% statistics
[h,p, t, df] =  ttestch(mean(L_energy), std(L_energy), L_energy(end), np-1, 0.05)
p

%% (3) load graphs, edge = midVOT
chlist = [1:25 28];
miVOT456 = load('mi_matrix_VOT456_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
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

diff_midVOT(diff_midVOT<0)=0; %< is the right one
diff_midVOT = abs(diff_midVOT);

energy_Stime = zeros(np,nwin);
energy_Xtime = zeros(np,nwin);
S_pat_m = zeros(np,nch,nch,nwin);
L_pat_m = zeros(np,nch,nch);
parfor ipat = 1:np
    tic
    X = squeeze(diff_midVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    S_pat_m(ipat,:,:,:)=Sx;
    L_pat_m(ipat,:,:)=Lx;
    toc
end



mid_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(S_pat_m(ipat,:,:,iwin));
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
%title('Energy of matrices: MidVOT condition');
set(gcf,'color','w');
box off;
ax=gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;


figure;
for ipat = 1:np
     subplot(5,5,ipat);
     imagesc(squeeze(L_pat_m(ipat,:,:)));
     caxis([0 0.03]);

end

figure('position',[440,378,900,420]);
%figure;
L_energy = zeros(1,np);
for ipat = 1:np
    temp = squeeze(L_pat_m(ipat,:,:));
    temp(temp<0)=0;
    L_energy(ipat) = norm(temp,'fro');
end
h = bar(1:21,L_energy(1:21),'b');
hold on;
h1 = bar(22,L_energy(22),'r');
set(gcf,'color','w');
box off;
ax=gca;
ax.XTick = 1:22;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

% statistics
[h,p, t, df] =  ttestch(mean(L_energy), std(L_energy), L_energy(end), np-1, 0.05)
p


%% statistics for low S matrices

X1 = low_VOT_energy;

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



%% statistics for high S matrices

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



%% statistics for mid S matrices

X1 = mid_VOT_energy;

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



