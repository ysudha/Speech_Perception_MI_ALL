% let's vary n, how do the results compare?
%% single graph
nlist = 5:50;
nWin = 5;
normMat = zeros(length(nlist),nWin);
for in = 1:length(nlist)
n = nlist(in);
p = 0.15; %0.09 works

rand('seed',180); % reseed so you get a similar picture %100, and 180 work
G = rand(n,n) < p;

origLowRank = G;
origLowRank = triu(origLowRank,1);
origLowRank = origLowRank + origLowRank';
% figure;
% imagesc(origLowRank);

%% use this as the core, and add another matrix on top, in 40 time windows

nWin = 5;
simMat = zeros(n,n,nWin);

pSparse = 0.12;
for i = 1:nWin
    Gtemp = (rand(n,n)<pSparse) + G;
    Gtemp(Gtemp>0) = 1;
    Gtemp = triu(Gtemp,1);
    Gtemp = Gtemp + Gtemp';
    simMat(:,:,i) = Gtemp;
end

%% let's plot these matrices
% 
% figure;
% for i = 1:nWin
%     subplot(2,5,i);
%     Gr = graph(squeeze(simMat(:,:,i)));
%     plot(Gr,'Layout','circle','NodeFontSize',12,'LineWidth',2,'MarkerSize',8,'EdgeColor','k','EdgeAlpha',0.6);
%     axis square;
%     box off;
%     axis off;
% end
% 
% 
% for i = 1:nWin
%     subplot(2,5,i+5);
%     imagesc(squeeze(simMat(:,:,i)));
%     axis square;
% end
% set(gcf,'color','w');

%% robust PCA time 
[L,S] = RobustPCA_time(simMat(:,:,:));
% figure;
% for i = 1:nWin
%     subplot(1,5,i);
%     imagesc(squeeze(simMat(:,:,i)));
%      axis square;
% end
% 
% figure;
% for i = 1:nWin
%     subplot(1,5,i);
%     t = squeeze(S(:,:,i));
%     t(t<0)=0;
%     imagesc(t);
%      axis square;
% end
% figure;
% imagesc(L);


for i = 1:nWin
   normMat(in,i) = norm((squeeze(simMat(:,:,i))-L-S(:,:,i)),'fro');
end

end


figure; plot(mean(normMat,2))
