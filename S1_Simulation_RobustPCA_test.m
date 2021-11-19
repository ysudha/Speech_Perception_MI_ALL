%% single graph
n = 8;
p = 0.09; %0.09 + 180 works. only 2 connections

rand('seed',180); % reseed so you get a similar picture %100, and 180 work
G = rand(n,n) < p;

origLowRank = G;
origLowRank = triu(origLowRank,1);
origLowRank = origLowRank + origLowRank';
figure;
imagesc(origLowRank);

%% use this as the core, and add another matrix on top, in 40 time windows

nWin = 5;
simMat = zeros(n,n,nWin);
rand('seed',180);
pSparse = 0.12;
for i = 1:nWin
    Gtemp = (rand(n,n)<pSparse) + G;
    Gtemp(Gtemp>0) = 1;
    Gtemp = triu(Gtemp,1);
    Gtemp = Gtemp + Gtemp';
    simMat(:,:,i) = Gtemp;
end

%% let's plot these matrices

figure;
for i = 1:nWin
    subplot(2,nWin,i);
    Gr = graph(squeeze(simMat(:,:,i)));
    plot(Gr,'Layout','circle','NodeFontSize',30,'LineWidth',2,'MarkerSize',8,'EdgeColor','k','EdgeAlpha',0.6);
    axis square;
    box off;
    axis off;
end

for i = 1:nWin
    subplot(2,nWin,i+nWin);
    imagesc(squeeze(simMat(:,:,i)));
    axis square;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',20,'FontWeight','bold')
    set(gca,'XTickLabelMode','auto')
end
set(gcf,'color','w');

%% robust PCA time 
[L,S] = RobustPCA_time(simMat(:,:,:));
figure;
for i = 1:nWin
    subplot(1,nWin,i);
    imagesc(squeeze(simMat(:,:,i)));
     axis square;
end

figure;
for i = 1:nWin
    subplot(2,nWin,i);
    t = squeeze(S(:,:,i));
    %t(t<0)=0;
    imagesc(t);
     axis square;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
set(gcf,'color','w');
end

figure;
for i = 1:nWin
    subplot(2,nWin,i);
    t = squeeze(S(:,:,i));
    t(t<0.1)=0;
    Gr = graph(t,'upper');
    plot(Gr,'Layout','circle','NodeFontSize',30,'LineWidth',2,'MarkerSize',8,'EdgeColor','k','EdgeAlpha',0.6);
    axis square;
    box off;
    axis off;
end
set(gcf,'color','w');


figure;
subplot(2,nWin,1);
imagesc(L);
L(L<0.01)=0;
axis square;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
set(gcf,'color','w');

figure;
for i = 1
    subplot(2,5,i);
    Gtemp = triu(L,1);
    Gtemp = L + L';
    Gr = graph(Gtemp);
    plot(Gr,'Layout','circle','NodeFontSize',30,'LineWidth',2,'MarkerSize',8,'EdgeColor','k','EdgeAlpha',0.6);
    axis square;
    box off;
    axis off;
end

normMat = zeros(1,nWin);
norm_orig = zeros(1,nWin);
for i = 1:nWin
   normMat(i) = norm((squeeze(simMat(:,:,i))-L-S(:,:,i)),'fro');
  norm_orig(i) = norm((squeeze(simMat(:,:,i))),'fro');
end
figure;
plot(normMat); hold on;
plot(norm_orig);

% %% running robust PCA with time component
% 
% %[L,S] = RobustPCA_time(simMat(:,:,1:end));
% temp = rand(n,n,nWin+1);
% for i = 1:nWin+1
%     Gtemp = squeeze(temp(:,:,i));
%    Gtemp = triu(Gtemp,1);
%     Gtemp = Gtemp + Gtemp';
%     temp(:,:,i) = Gtemp;
% end
% for i = 2:nWin+1
%    temp(:,:,i-1) = temp(:,:,i) -temp(:,:,1);
% end
% 
% temp = temp(:,:,2:nWin+1);
% 
% [L,S] = RobustPCA_time(temp);
% figure;
% for i = 1:nWin
%     subplot(8,5,i);
%     imagesc(squeeze(temp(:,:,i)));
% end
% 
% figure;
% for i = 1:nWin
%     subplot(8,5,i);
%     imagesc(squeeze(S(:,:,i)));
% end
% figure;
% imagesc(L);
% 
% normMat = zeros(1,nWin);
% norm_orig = zeros(1,nWin);
% for i = 1:nWin
%    normMat(i) = norm((squeeze(temp(:,:,i))-L-S(:,:,i)),'fro');
%   norm_orig(i) = norm((squeeze(temp(:,:,i))),'fro');
% end
% figure;
% plot(normMat); hold on;
% plot(norm_orig);
% 
% %% given the original G, how to plot
% 
