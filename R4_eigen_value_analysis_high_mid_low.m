% (1) load graphs, edge = midVOT - lowVOT ...looking at increases and decreases
% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases

k=10; % lets assume rank k is fixed for now

% excluding pat#10 improves PCA... very susceptible to outliers after all
 %% (1) load graphs, edge = midVOT - lowVOT ...looking at increases and decreases
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
error_curves = zeros(np,nwin);
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
%     figure;
%     plot(3:nwin,error_curves(ipat,3:end),'-*','LineWidth',4);
    xlabel('Time');
    ylabel('Approximation error, using rank k apprx');
    set(gcf,'color','w'); box off;
end

% plotting
X1 = error_curves(:,3:end);
X1(10,:)=[];
[coeff, score, latent, tsquared, explained, mu] = pca(X1(:,1:end));%20
figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');

Xhat1 = score(:,1:3)*coeff(:,1:3)' + mu;
 figure; imagesc(X1); figure; imagesc(Xhat1)
 
 plot(Xhat1(1:end-1,:)','-b');hold on;
plot(Xhat1(21,:)','-r');hold on;

%% (2) load graphs, edge = midVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = midVOT - highVOT ...looking at increases only
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin);
error_curves = zeros(np,nwin);
for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 1:nwin
       mat = squeeze(squeeze(mi456{iwin,1} - mi789{iwin,1}) );
       %mat(mat<0)=0;
       diff_mid_lowVOT(ipat,:,:,iwin)= mat;
       [U,S,V] = eig(mat);
       apprx_err = norm(mat - U(:,1:k)*S(1:k,1:k)*V(:,1:k)','fro');
       error_curves(ipat,iwin) = apprx_err;
    end
%     figure;
%     plot(3:nwin,error_curves(ipat,3:end),'-*','LineWidth',4);
    xlabel('Time');
    ylabel('Approximation error, using rank k apprx');
    set(gcf,'color','w'); box off;
end

% plotting
X2 = error_curves(:,3:end);
X2(10,:)=[];
[coeff, score, latent, tsquared, explained, mu] = pca(X2(:,1:end));%20
figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');

Xhat2 = score(:,1:3)*coeff(:,1:3)' + mu;
 figure; imagesc(X2); figure; imagesc(Xhat2)
 
 plot(Xhat2(1:end-1,:)','-b');hold on;
plot(Xhat2(21,:)','-r');hold on;

%% (3) load graphs, edge = lowVOT - highVOT ...looking at increases and decreases
%%%% load new graphs, edge = lowVOT - highVOT ...looking at increases only
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin);
error_curves = zeros(np,nwin);
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 1:nwin
       mat = squeeze(squeeze(mi456{iwin,1} - mi123{iwin,1}) );
       %mat(mat<0)=0;
       diff_mid_lowVOT(ipat,:,:,iwin)= mat;
       [U,S,V] = eig(mat);
       apprx_err = norm(mat - U(:,1:k)*S(1:k,1:k)*V(:,1:k)','fro');
       error_curves(ipat,iwin) = apprx_err;
    end
%     figure;
%     plot(3:nwin,error_curves(ipat,3:end),'-*','LineWidth',4);
    xlabel('Time');
    ylabel('Approximation error, using rank k apprx');
    set(gcf,'color','w'); box off;
end

% plotting
X3 = error_curves(:,3:end);
X3(10,:)=[];
[coeff, score, latent, tsquared, explained, mu] = pca(X3(:,1:end));%20
figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');

Xhat3 = score(:,1:3)*coeff(:,1:3)' + mu;
 figure; imagesc(X3); figure; imagesc(Xhat3)
 
 plot(Xhat3(1:end-1,:)','-b');hold on;
plot(Xhat3(21,:)','-r');hold on;

%% Combine the error vectors from all conditions

X = [X1 X2 X3];

[coeff, score, latent, tsquared, explained, mu] = pca(X(:,1:end));%20

figure; scatter(score(1:end-1,1),score(1:end-1,2),'b*');hold on;scatter(score(end,1),score(end,2),'r*');
figure; scatter3(score(1:end-1,1),score(1:end-1,2),score(1:end-1,3),'b*');hold on;scatter3(score(end,1),score(end,2),score(end,3),'r*');

Xhat = score(:,1:3)*coeff(:,1:3)' + mu;
 figure; imagesc(X); figure; imagesc(Xhat)
 
 plot(Xhat(1:end-1,:)','-b');hold on;
plot(Xhat(21,:)','-r');hold on;

