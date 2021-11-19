% let's vary number of windows, how do the results compare?
%% single graph
nlist =5:75;
n=8; % fixed n
normMat = zeros(1,length(nlist));
normOrig = zeros(1,length(nlist));
std_normMat = zeros(1,length(nlist));
std_normOrig = zeros(1,length(nlist));
for in = 1:length(nlist)
    n = 8;
    p = 0.15; %0.09 works

    rand('seed',180); % reseed so you get a similar picture %100, and 180 work
    G = rand(n,n) < p;

    origLowRank = G;
    origLowRank = triu(origLowRank,1);
    origLowRank = origLowRank + origLowRank';
    % figure;
    % imagesc(origLowRank);

    %% use this as the core, and add another matrix on top, in 40 time windows

    nWin = nlist(in);
    simMat = zeros(n,n,nWin);

    pSparse = 0.12;
    for i = 1:nWin
        Gtemp = (rand(n,n)<pSparse) + G;
        Gtemp(Gtemp>0) = 1;
        Gtemp = triu(Gtemp,1);
        Gtemp = Gtemp + Gtemp';
        simMat(:,:,i) = Gtemp;
    end



    %% robust PCA time 
    [L,S] = RobustPCA_time(simMat(:,:,:));

    % for i = 1:nWin
    %     t = squeeze(S(:,:,i));
    %     t(t<0)=0;
    %     S(:,:,i)=t;
    % end

    tempnormMat = zeros(1,nWin);
    norm_orig = zeros(1,nWin);
    for i = 1:nWin
       tempnormMat(i) = (norm((squeeze(simMat(:,:,i))-L-S(:,:,i)),'fro'));
         norm_orig(i) = norm((squeeze(simMat(:,:,i))),'fro');
    end

    normMat(in)=  mean(tempnormMat);
    normOrig(in)=  mean(norm_orig);

    std_normMat(in)=  tinv(0.975,nWin-1)*std(tempnormMat)./sqrt(nWin);
    std_normOrig(in)=  tinv(0.975,nWin-1)*std(norm_orig)./sqrt(nWin);
end


%%
figure; errorbar(nlist(10:end),normMat(10:end),std_normMat(10:end)); hold on;
ylim([0,0.001]);
%errorbar(nlist,(normOrig),std_normOrig );
