% (1) plots edges as increase of mid- low
% (2) plots edges as decrease of mid-low
% (3) plots edges as mid-low, considering both + and _ values, ignoring
% person 10, outlier was the best bet
%% (1) load new graphs, edge = midVOT - lowVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin);
error_curves = zeros(ipat,nwin);
k=10; % lets assume rank k is fixed for now
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 1:nwin
       mat = squeeze(mi456{iwin,1}) - squeeze(mi123{iwin,1});
       mat(mat<0)=0;
       diff_mid_lowVOT(ipat,:,:,iwin)= mat;
       [U,S,V] = eig(mat);
       apprx_err = norm(mat - U(:,1:k)*S(1:k,1:k)*V(:,1:k)','fro');
       error_curves(ipat,iwin) = apprx_err;
    end
    figure;
    plot(3:nwin,error_curves(ipat,3:end),'-*','LineWidth',4);
    xlabel('Time');
    ylabel('Approximation error, using rank k apprx');
    set(gcf,'color','w'); box off;
end

%% (2)load graphs, edge = midVOT - lowVOT ...looking at decreases only
X = error_curves(:,3:end);
[coeff, score, latent, tsquared, explained, mu] = pca(X(:,1:20));%20
figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
%figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter(score(end,1),score(end,2),score(end,3),'r*');

Xhat = score(:,1)*coeff(:,1)' + mu;

%%%% load new graphs, edge = midVOT - lowVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin);
error_curves = zeros(ipat,nwin);
k=10; % lets assume rank k is fixed for now
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 1:nwin
       mat = squeeze(squeeze(mi123{iwin,1} - mi456{iwin,1}) );
       mat(mat<0)=0;
       diff_mid_lowVOT(ipat,:,:,iwin)= mat;
       [U,S,V] = eig(mat);
       apprx_err = norm(mat - U(:,1:k)*S(1:k,1:k)*V(:,1:k)','fro');
       error_curves(ipat,iwin) = apprx_err;
    end
    figure;
    plot(3:nwin,error_curves(ipat,3:end),'-*','LineWidth',4);
    xlabel('Time');
    ylabel('Approximation error, using rank k apprx');
    set(gcf,'color','w'); box off;
end

%%
X = error_curves(:,3:end);
[coeff, score, latent, tsquared, explained, mu] = pca(X(:,1:20));%20
figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
%figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter(score(end,1),score(end,2),score(end,3),'r*');

Xhat = score(:,1)*coeff(:,1)' + mu;
 figure; imagesc(X); figure; imagesc(Xhat)
 
 
 %% (3) load graphs, edge = midVOT - lowVOT ...looking at increases and decreases
X = error_curves(:,3:end);
[coeff, score, latent, tsquared, explained, mu] = pca(X(:,1:20));%20
figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
%figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter(score(end,1),score(end,2),score(end,3),'r*');

Xhat = score(:,1)*coeff(:,1)' + mu;

%%%% load new graphs, edge = midVOT - lowVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin);
error_curves = zeros(ipat,nwin);
k=10; % lets assume rank k is fixed for now
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 1:nwin
       mat = squeeze(squeeze(mi456{iwin,1} - mi123{iwin,1}) );
 %      mat(mat<0)=0;
       diff_mid_lowVOT(ipat,:,:,iwin)= mat;
       [U,S,V] = eig(mat);
       apprx_err = norm(mat - U(:,1:k)*S(1:k,1:k)*V(:,1:k)','fro');
       error_curves(ipat,iwin) = apprx_err;
    end
    figure;
    plot(3:nwin,error_curves(ipat,3:end),'-*','LineWidth',4);
    xlabel('Time');
    ylabel('Approximation error, using rank k apprx');
    set(gcf,'color','w'); box off;
end

%%
X = error_curves(:,3:end);
X(10,:)=[];
[coeff, score, latent, tsquared, explained, mu] = pca(X(:,1:end));%20
figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');

Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
 figure; imagesc(X); figure; imagesc(Xhat)
 
 plot(Xhat(1:end-1,:)','-b');hold on;
plot(Xhat(21,:)','-r');hold on;