% this code is similar to R19. Looking at which regions of the brain make
% Low-mid >0 and High-mid>0 the bext edge metrics


chlist = [1:25 28];%[1 3 4 8 12 13 2 6 7 11 15 16];%[1:25 28];
 %% (1) load graphs, edge = lowVOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Robust PCA in time%%%%%%%%%%%%%%%%%%

miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_lowVOT = zeros(np,nch,nch,nwin-2);
diff_mid_lowVOT = zeros(np,nch,nch,nwin-2);

for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 3:nwin
       mat1 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       mat = squeeze(mat1-mat2);
       diff_mid_lowVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

S_mat = cell(1,np);
parfor ipat = 1:np
    tic
    Y = squeeze(diff_mid_lowVOT(ipat,:,:,:));
    Y(Y>0)=0;
    Y=Y.*(-1);
    [Ly,Sy] = RobustPCA_time(Y);
    S_mat{1,ipat} = Sy;
    toc
end

% plots to look at frobenius borm of S in all windows.. comparing patients.
% Some windows clearly have greater energy in the matrix for aphasia
% patient
for i = 1:41
    figure;hold on;
    
    for ipat = 1:np
        
        S = S_mat{1,ipat};
        %imagesc(squeeze(S(:,:,ipat)));
        %hist(sum(squeeze(S(:,:,ipat))))
        bar(ipat,norm(squeeze(S(:,:,i)),'fro'));
    end
end


diff_mid_lowVOT(diff_mid_lowVOT>0)=0;
diff_mid_lowVOT = abs(diff_mid_lowVOT);

for ipat = 1:np
    figure;
    mat = squeeze(diff_mid_lowVOT(ipat,:,:,1));
    subplot(1,2,1);imagesc(mat);
    %subplot(1,2,2);hist(sum(mat));
    %C=clustering_coef_wu(mat)
    subplot(1,2,2);%imagesc(C);
     bar(ipat,norm(mat,'fro'));
    %str = sprintf('%0.2f',mean(C));
    %title(str);
    norm(mat,'fro')
    
end


for i = 1:41
    figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_mid_lowVOT(ipat,:,:,i));
        bar(ipat,norm(mat,'fro'));
       
    end
end

% looking at magnitudes with threshold
thresh = 0.2;
for ipat = 1:np
    figure;
%     S = S_mat{1,ipat};
%     S = squeeze(S(:,:,1));
    mat = squeeze(diff_mid_lowVOT(ipat,:,:,1));
    mat(mat<thresh)=0;
    G = graph(mat);
%     S(S<thresh)=0;
%     G = digraph(S);
    plot(G);
end


% looking at degree distribution %does not differentiate the nodes

%% edges in S with the greatest strength
% find the edges with greatest strength
mask_all = zeros(nch,nch,nwin);
for iwin = 1:41
    figure;
    mask = ones(nch,nch);
    for ipat = 1:np-1
        S = S_mat{1,ipat};
        S = squeeze(S(:,:,iwin));
        Sap = S_mat{1,22};
        Sap = squeeze(Sap(:,:,iwin));
        mask_subj_specific = Sap>S;
        mask = mask.*mask_subj_specific;
    end
    imagesc(mask)
    mask_all(:,:,iwin) = mask;
end

% relevant time windows which were significant
sig_win = [1:4 23:31 33 37:41];
sig_win = [1:4];

mask = squeeze(mask_all(:,:,sig_win(1)));
for i = 2:length(sig_win)
    mask = mask.*squeeze(mask_all(:,:,sig_win(i)));
end
figure;
imagesc(mask)

%% edges in the actual matrix
% find the edges with greatest strength
mask_all = zeros(nch,nch,nwin);
for iwin = 1:41
    figure; hold on;
    mask = ones(nch,nch);
    for ipat = 1:np-1
        S = squeeze(diff_mid_lowVOT(ipat,:,:,iwin));
        Sap = squeeze(diff_mid_lowVOT(22,:,:,iwin));
        mask_subj_specific = Sap>S;
        mask = mask.*mask_subj_specific;
    end
    imagesc(mask)
    mask_all(:,:,iwin) = mask;
end

% relevant time windows which were significant
sig_win = [1:4 23:31 33 37:41];
sig_win = [1:4];

mask = squeeze(mask_all(:,:,sig_win(1)));
for i = 2:length(sig_win)
    mask = mask.*squeeze(mask_all(:,:,sig_win(i)));
end
figure;
imagesc(mask)





