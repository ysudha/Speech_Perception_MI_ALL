% combining areas to see how well the aphasia patient separates 
% using best edge matrices, mid-low, high-low, mid-high, baseline
% normalized
% (1) load graphs, edge = lowVOT ...looking at increases and decreases
% (2) load graphs, edge = midVOT ...looking at increases and decreases
% (3) load graphs, edge = highVOT ...looking at increases and decreases
%% regions to be grouped
% [1,2] [ 3,4], [5,9,10,14], [6,7],
% [8,13,18],[11,15,19],[12],[16],[17,21],[20,25], [22,23,24,26,27,28,29,30]
grp = {[1,2], [ 3,4], [5,9,10,14], [6,7],[8,13,18],[11,15,19],[12],[16],[17,21],[20,25],[22,23,24,28]};%, [22,23,24,26,27,28,29,30]};

%% (1) load graphs, edge = lowVOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%

miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(grp);
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin-2);
diff_low_midVOT = zeros(np,nch,nch,nwin-2);
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       mat_m_l = squeeze(mat1-mat2);
       mat_m_l(mat_m_l<0)=0;
       mat_l_m = squeeze(mat2-mat1);
       mat_l_m(mat_l_m<0)=0;
       diff_mid_lowVOT(ipat,:,:,iwin-2) = contrct(mat_m_l,grp);
       diff_low_midVOT(ipat,:,:,iwin-2) = contrct(mat_l_m,grp);
    end

end

nwin = nwin-2;
energy_S_ml_time = zeros(np,nwin);
energy_S_lm_time = zeros(np,nwin);
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_low_midVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_ml_time(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_lm_time(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_ml_time(:,1:end);
X2 = energy_S_lm_time(:,1:end);

figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-low >0');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-low <0');


[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid>low');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid<low');


%% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(grp);
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT789, only
% positive increases diff_mid_highVOT
diff_mid_highVOT = zeros(np,nch,nch,nwin-2);
diff_high_midVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} -mi456{1,1} ); 
       mat2 = squeeze(mi789{iwin,1}- mi789{1,1}); 
       mat_m_h = squeeze(mat1-mat2);
       mat_m_h(mat_m_h<0)=0;
       mat_h_m = squeeze(mat2-mat1);
       mat_h_m(mat_h_m<0)=0;
       diff_mid_highVOT(ipat,:,:,iwin-2) = contrct(mat_m_h,grp);
       diff_high_midVOT(ipat,:,:,iwin-2) = contrct(mat_h_m,grp);
    end

end

nwin = nwin-2;
energy_S_mh_time = zeros(np,nwin);
energy_S_hm_time = zeros(np,nwin);
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_high_midVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mh_time(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_hm_time(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_mh_time(:,1:end);
X2 = energy_S_hm_time(:,1:end);

figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-high >0');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-high <0');


[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid>high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid<high');


%% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT789 = load('mi_matrix_VOT789_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(grp);
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT789, only
% positive increases diff_low_highVOT
diff_low_highVOT = zeros(np,nch,nch,nwin-2);
diff_high_lowVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi123{iwin,1}- mi123{1,1}); 
       mat2 = squeeze(mi789{iwin,1}- mi789{1,1}); 
       mat_l_h = squeeze(mat1-mat2);
       mat_l_h(mat_l_h<0)=0;
       mat_h_l = squeeze(mat2-mat1);
       mat_h_l(mat_h_l<0)=0;
       diff_low_highVOT(ipat,:,:,iwin-2) = contrct(mat_l_h,grp);
       diff_high_lowVOT(ipat,:,:,iwin-2) = contrct(mat_h_l,grp);   
    end

end

nwin = nwin-2;
energy_S_lh_time = zeros(np,nwin);
energy_S_hl_time = zeros(np,nwin);
parfor ipat = 1:np
    tic
    X = squeeze(diff_low_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_high_lowVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_lh_time(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_hl_time(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_lh_time(:,1:end);
X2 = energy_S_hl_time(:,1:end);

figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high >0');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high <0');


[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low>high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low<high');

%% plotting 
figure; hold on;
for ipat = 1:np

    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    for iwin = 20%1:nwin
        subplot(5,5,ipat);
        mat = squeeze(X(:,:,iwin));
        imagesc(mat);
        colormap jet; caxis([0 0.3]);
        axis square;
        str = sprintf('%0.2f',energy_S_ml_time(ipat,iwin));
        title(str);
    end
end

figure; hold on;
for ipat = 1:np

    X = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(X);
    for iwin = 20%1:nwin
        subplot(5,5,ipat);
        mat = squeeze(Sy(:,:,iwin));
        [U1,S1,V1] = eig(mat);
        imagesc(U1(:,1:1)*S1(1:1,1:1)*V1(:,1:1)');
        colormap jet; caxis([0 0.7]);
        axis square;
        str = sprintf('%0.2f',energy_S_ml_time(ipat,iwin));
        title(str);
    end
end
%%  function to contract matrix into regions
function [regmat] = contrct(mat, grp)
    ngrp = length(grp);
    regmat = zeros(ngrp,ngrp);
    for ir = 1:ngrp
        for ic = 1:ngrp
            vecr = grp{ir}; %list of row indices for ir row in region matrix
            vecc = grp{ic}; %list of col indices for ic col in region matrix
            submatrix = mat(vecr,vecc);
            if ir==ic
                regmat(ir,ic) = sum(diag(submatrix)) + sum(sum(submatrix - diag(diag(submatrix))))/2; %the submatrix is symmetric. Keep the diagonals separate, add to the just the mean of the sum of the remaining elements. (One of the triangular elements)
            else
                regmat(ir,ic) = sum(sum(submatrix)); %the submatrix is not symmetric. add all the elements.
            end
        end
    end
end %end function



