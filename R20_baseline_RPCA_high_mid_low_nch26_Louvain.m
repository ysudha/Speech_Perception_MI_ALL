% same as R18, with  fewer channels, excluding 26,27,29,30, community
% detection

chlist = [1:25 28];%1:30; %[1 3 4 8 12 13 ];%[2 6 7 11 15 16];%[1:25 28];
%% (1) load graphs, edge = BSmid - BSlowVOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%

miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
n=nch;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       mat =  squeeze(mat1-mat2);
       diff_mid_lowVOT(ipat,:,:,iwin-2) =mat(chlist, chlist);
    end
end

nwin = 41;

%%%%%%%%%normalization, keeping only positive values
% for ipat = 1:np
%     for iwin = 1:nwin
%             mat = squeeze(diff_mid_lowVOT(ipat,:,:,iwin));
%             mat(mat<0)=0;
%             diff_mid_lowVOT(ipat,:,:,iwin) = mat;
%     end
% end

rnd_diff_mid_lowVOT = diff_mid_lowVOT;
Qmatrix = zeros(np,nwin);
for ipat = 5
    mat = squeeze(diff_mid_lowVOT(ipat,:,:,:));
%random matrix, and recalculate the clustering coefficient.
for iwin = 1%1:nwin
  %  mat = squeeze(diff_mid_lowVOT(ipat,:,:,iwin));
%           C = clustering_coef_wu(mat);
% 
%         uppertriangle= find(triu(mat,1)>0);
%         matrixr = zeros(size(mat));
%         matrixr(uppertriangle) = mat(uppertriangle(randperm(numel(uppertriangle)))); 
%         matrixr = matrixr+matrixr';
%         Cr = clustering_coef_wu(matrixr);
%         subplot(7,7,iwin); bar(mean(C'),mean(Cr'));

        % Let's look at the distribution of clustering coefficients in the random
        % network overlad on the distribution from the real network.

%         subplot(4,4,iwin-14); hist([C'; Cr']'); 
%         str = sprintf('%0.3f, %0.3f',mean(C),mean(Cr));
%         title(str);
        [L,S] = RobustPCA_time(mat);
        matrixr = zeros(nch,nch);
        [partition1,Q1]=community_louvain(squeeze(mat(:,:,iwin)),1,randperm(n),'negative_asym');
        
        [Uq,Sq,Vq] = eig(squeeze(S(:,:,iwin)));
        energy_S_mat = sum(abs(diag(Sq)))
        [indSort,CiSort] = fcn_order_partition(mat,partition1);
        [x,y]=fcn_grid_communities(CiSort);
        figure; imagesc(mat(indSort,indSort)); hold on; plot(x,y,'k');axis square;
        
        

        Qmatrix(ipat,iwin) = Q1;
        nperm = 2;
        Qran = zeros(1,nperm);
        Eran = zeros(1,nperm);
        for j = 1:nperm
            matrand = zeros(size(mat));
            for ii = 1:nwin
                m = squeeze(mat(:,:,ii));
                matrixr = zeros(n,n);
                uppertriangle= find(triu(m,1)~=0);
                matrixr(uppertriangle) = m(uppertriangle(randperm(numel(uppertriangle))));
                matrixr = matrixr+matrixr';
                matrand(:,:,ii) = matrixr;
            end
            [L,S] = RobustPCA_time(matrand);
            mr = squeeze(matrand(:,:,iwin));

            [partitionr1,Qr1]=community_louvain(mr,1,randperm(n),'negative_asym');
            Qran(j)  = Qr1;
            [Uq,Sq,Vq] = eig(squeeze(S(:,:,iwin)));
            Eran(j) = sum(abs(diag(Sq)));
            [indSort,CiSort] = fcn_order_partition(squeeze(matrand(:,:,iwin)),partitionr1);
            [x,y]=fcn_grid_communities(CiSort);
            figure; imagesc(mr(indSort,indSort)); hold on; plot(x,y,'k');axis square;
            
        
            
        end
%         subplot(7,7,iwin)
%         errorbar([1],[mean(Qran)],[std(Qran)],'s','MarkerFaceColor','red');
%         hold on; plot(2,Q1,'b*');
%         figure;
%         [indSort,CiSort] = fcn_order_partition(mat,partition1);
%         [x,y]=fcn_grid_communities(CiSort);
%         subplot(1,2,1); imagesc(mat(indSort,indSort)); hold on; plot(x,y,'k');
% 
%         %And we can compare to what we would see in the random graph like this:
% 
%         [indSort,CiSort] = fcn_order_partition(matrixr,partitionr1);
%         [x,y]=fcn_grid_communities(CiSort);
%         subplot(1,2,2); imagesc(matrixr(indSort,indSort)); hold on; plot(x,y,'k');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %randomize matrices, how ddoes it function
% for iwin = 1:nwin
%     mat = squeeze(diff_mid_lowVOT(ipat,:,:,iwin));
%         uppertriangle= find(triu(mat,1)>0);
%         matrixr = zeros(size(mat));
%         matrixr(uppertriangle) = mat(uppertriangle(randperm(numel(uppertriangle)))); 
%         matrixr = matrixr+matrixr';
%         rnd_diff_mid_lowVOT(ipat,:,:,iwin) = matrixr;
% end
% mat = squeeze(rnd_diff_mid_lowVOT(ipat,:,:,:));
% % end randomization

% M.matrices = mat;
% gamma = 1.01;
% omega = 0.2;
% T=length(squeeze(M.matrices(1,1,:)));
% B=spalloc(n*T,n*T,n*n*T+2*n*T);
% twomu=0;
% for s=1:T
%     k=sum(squeeze(M.matrices(:,:,s)));
%     twom=sum(k);
%     twomu=twomu+twom;
%     indx=[1:n]+(s-1)*n;
%     B(indx,indx)=squeeze(M.matrices(:,:,s))-gamma*k'*k/twom;
% end
% twomu=twomu+2*omega*n*(T-1);
% B = B + omega*spdiags(ones(n*T,2),[-n,n],n*T,n*T);
% [S,Q] = genlouvain(B,10000,0);
% Q = Q/twomu;
% partition_multi = reshape(S,n,T);
% 
% % What does the partition of nodes into communities look like?
% 
% figure; imagesc(partition_multi(:,:))
% 
% Q
%     

end

figure; imagesc(Qmatrix)
    
%% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 26;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT789, only
% positive increases diff_mid_highVOT
diff_midVOT = zeros(np,nch,nch,nwin-2);
diff_mid_highVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} -mi456{1,1} ); 
       mat2 = squeeze(mi789{iwin,1}- mi789{1,1}); 
       diff_midVOT(ipat,:,:,iwin-2)= squeeze(mat1(chlist,chlist));
       mat = squeeze(mat1-mat2);
       diff_mid_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;
energy_S_mid_high = zeros(np,nwin);
energy_S_mid = zeros(np,nwin);
S_pat_mh_hm = cell(1,np);
L_pat_mh_hm = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_mid_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_midVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    S_pat_mh_hm{1,ipat}=Sx;
    L_pat_mh_hm{1,ipat}=Lx;
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_high(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_mid(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_mid_high(:,1:end);
X2 = energy_S_mid(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S mid-high');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S just mid');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid-high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = mid');

% reconstruct with 1st PC
[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
 figure;
 plot(Xhat(1:end-1,:)','-b');hold on;
plot(Xhat(end,:)','-r');hold on;

%% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT789 = load('mi_matrix_VOT789_K4.mat'); 
np = 22; %number of people end undergrads, 1 aphasia
nwin = 43;
nch = 26;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT789, only
% positive increases diff_low_highVOT
diff_highVOT = zeros(np,nch,nch,nwin-2);
diff_low_highVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi123{iwin,1}- mi123{1,1}); 
       mat2 = squeeze(mi789{iwin,1}- mi789{1,1}); 
       diff_highVOT(ipat,:,:,iwin-2)= squeeze(mat2(chlist,chlist));
       mat = squeeze(mat1-mat2);
       diff_low_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;
energy_S_low_high = zeros(np,nwin);
energy_S_high = zeros(np,nwin);
S_pat_lh_hl = cell(1,np);
L_pat_lh_hl = cell(1,np);
parfor ipat = 1:np
    tic
    X = squeeze(diff_low_highVOT(ipat,:,:,:));
    [Lx,Sx] = RobustPCA_time(X);
    Y = squeeze(diff_highVOT(ipat,:,:,:));
    [Ly,Sy] = RobustPCA_time(Y);
    S_pat_lh_hl{1,ipat}=Sx;
    L_pat_lh_hl{1,ipat}=Lx;
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_low_high(ipat,iwin) = sum(abs(diag(S1)));
        [U2,S2,V2] = eig(squeeze(Sy(:,:,iwin)));
        energy_S_high(ipat,iwin) = sum(abs(diag(S2)));
    end
    toc
end


X1 = energy_S_low_high(:,1:end);
X2 = energy_S_high(:,1:end);
figure; 
 plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high');

figure; 
 plot(X2(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X2(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S just high');

[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = low-high');

[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('baseline norm Edge = high');

% reconstruct with 1st PC
[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
 figure;
 plot(Xhat(1:end-1,:)','-b');hold on;
plot(Xhat(end,:)','-r');hold on;

%% save the L and S matrices for all 3 conditions.
BS_Rpca.S_pat_ml_lm = S_pat_ml_lm;
BS_Rpca.L_pat_ml_lm = L_pat_ml_lm;
BS_Rpca.S_pat_mh_hm = S_pat_mh_hm;
BS_Rpca.L_pat_mh_hm = L_pat_mh_hm;
BS_Rpca.S_pat_lh_hl = S_pat_lh_hl;
BS_Rpca.L_pat_lh_hl = L_pat_lh_hl;
    
save('BS_Rpca','BS_Rpca','-v7.3');
%% this program compares baseline normalized matrices using edges mid-low, etc.
load('BS_Rpca');
energy_S_low_high = zeros(np,nwin);
energy_L_low_high = zeros(np,1);
for ipat = 1:np
    Sx = BS_Rpca.S_pat_lh_hl{1,ipat};
    Lx = BS_Rpca.L_pat_lh_hl{1,ipat};
    [U2,S2,V2] = eig(squeeze(Lx(:,:)));
    energy_L_low_high(ipat,1) = sum(abs(diag(S2)));
    
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_low_high(ipat,iwin) = sum(abs(diag(S1)));
    end  
end

X1 = energy_S_low_high(:,1:end);
ntim = [101:5:301];

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high');

figure(1); plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',1);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p] = ttest(X1(pat_list,iwin),X1(ipat,iwin));
        stats(iwin) = p;
        
    end
    %stats = fdr(stats);
    fprintf('Patient number %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end

%%%%%%%%%%%%%%
energy_S_mid_low = zeros(np,nwin);
energy_L_mid_low = zeros(np,1);
for ipat = 1:np
    Sx = BS_Rpca.S_pat_ml_lm{1,ipat};
    Lx = BS_Rpca.L_pat_ml_lm{1,ipat};
    [U2,S2,V2] = eig(squeeze(Lx(:,:)));
    energy_L_mid_low(ipat,1) = sum(abs(diag(S2)));
    
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_low(ipat,iwin) = sum(abs(diag(S1)));
    end  
end

X1 = energy_S_mid_low(:,1:end);
ntim = [101:5:301];

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high');

figure(1); plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',1);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p] = ttest(X1(pat_list,iwin),X1(ipat,iwin));
        stats(iwin) = p;
        
    end
    %stats = fdr(stats);
    fprintf('Patient number %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
energy_S_mid_high = zeros(np,nwin);
energy_L_mid_high = zeros(np,1);
for ipat = 1:np
    Sx = BS_Rpca.S_pat_mh_hm{1,ipat};
    Lx = BS_Rpca.L_pat_mh_hm{1,ipat};
    [U2,S2,V2] = eig(squeeze(Lx(:,:)));
    energy_L_mid_high(ipat,1) = sum(abs(diag(S2)));
    
    for iwin = 1:nwin
        [U1,S1,V1] = eig(squeeze(Sx(:,:,iwin)));
        energy_S_mid_high(ipat,iwin) = sum(abs(diag(S1)));
    end  
end

X1 = energy_S_mid_high(:,1:end);
ntim = [101:5:301];

figure; 
plot(X1(1:end-1,:)','-b*','LineWidth',4);hold on;
plot(X1(end,:)','-r*','LineWidth',4);hold on;
title('baseline norm Energy of S low-high');

figure(1); plot(ntim,mean(X1(1:end,:)),'*-k','LineWidth',1);hold on;
plot(ntim, X1(end,:),'-r*','LineWidth',1);
[ph,msg]=jbfill(ntim,(mean(X1(1:end,:))+std(X1(1:end,:))),(mean(X1(1:end,:))-std(X1(1:end,:))),'b','b',1,0.3);

stats = zeros(1,nwin);

for ipat = 1:np
    pat_list = 1:np;
    pat_list(ipat)=[];
    for iwin = 1:nwin
        [h,p] = ttest(X1(pat_list,iwin),X1(ipat,iwin));
        stats(iwin) = p;
        
    end
    stats = fdr(stats);
    fprintf('Patient number %d, numb of p<0.05 = %d\n',ipat, sum(stats<0.05));

end
    
