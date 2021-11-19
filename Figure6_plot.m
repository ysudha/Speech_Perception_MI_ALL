% time varying graph clustering for aphasic - vs healthy for low, mid and
% high VOT
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

% remove the negative values after bias correction
for ipat = 1:np
    for iwin = 1:nwin
        mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1};
        mi123(mi123<0)=0;
        miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1} = mi123;
    end
end

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
    if i~=1
        a=plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    else
        plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
    end
    hold on;
end
plot(ntim,low_VOT_energy(np,:),'r-*','LineWidth',4);


legend('Healthy subject','Aphasic subject');
ylabel('Frobenius norm');
xlabel('Time');
title('Energy of matrices: (LowVOT- BS) >0');
       
axis square;
set(gcf,'color','w');
box off;

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

%% combined LOWVOT figure
%let's plot the energy of baseline normalized graphs
figure('Position',[100 100 750 500]);
for i = 1:np-1
    if i~=1
        a=plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    else
        plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
    end
    hold on;
end
plot(ntim,low_VOT_energy(np,:),'r-o','LineWidth',4);
set(gca,'FontSize',16);
xlim([-25,425]);

%legend('Healthy subject','Aphasic subject','Location','bestoutside');
%legend box off;

ylabel('Frobenius norm','FontSize',20);
xlabel('Time (ms)','FontSize',20);

%title('Energy of matrices: (LowVOT- BS) >0');      

set(gcf,'color','w');
box off;
yl = ylim;

for i = 1:length(ntim)
    patstats = stats_low_gt_bs(22,:);
    if(patstats(i)==1)
        a=patch([ntim(i)-0.5*10 ntim(i)+0.5*10 ntim(i)+0.5*10 ntim(i)-0.5*10], [yl(2)-0.1 yl(2)-0.1 yl(2)-0.05 yl(2)-0.05],'y')
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';

        alpha(0.3);
        str = sprintf('*');
        text(ntim(i)-3,yl(2)-0.03,str,'FontSize',20);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
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

% remove the negative values after bias correction
for ipat = 1:np
    for iwin = 1:nwin
        mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1};
        mi789(mi789<0)=0;
        miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1} = mi789;
    end
end

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

high_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_highVOT(ipat,:,:,iwin));
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
title('Energy of matrices: (HighVOT- BS) >0');
%title('Energy of matrices: (LowVOT- BS) <0');

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

%% combined HighVOT figure
%let's plot the energy of baseline normalized graphs
figure('Position',[100 100 750 500]);
for i = 1:np-1
    if i~=1
        a=plot(ntim,high_VOT_energy(i,:),'b-*','LineWidth',4);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    else
        plot(ntim,high_VOT_energy(i,:),'b-*','LineWidth',4);
    end
    hold on;
end
plot(ntim,high_VOT_energy(np,:),'r-*','LineWidth',4);
set(gca,'FontSize',16);
xlim([-25,425]);

%legend('Healthy subject','Aphasic subject','Location','bestoutside');
%legend box off;

ylabel('Frobenius norm','FontSize',20);
xlabel('Time (ms)','FontSize',20);

%title('Energy of matrices: (HighVOT- BS) >0');      

set(gcf,'color','w');
box off;
yl = ylim;

for i = 1:length(ntim)
    patstats = stats_high_gt_bs(22,:);
    if(patstats(i)==1)
        a=patch([ntim(i)-0.5*10 ntim(i)+0.5*10 ntim(i)+0.5*10 ntim(i)-0.5*10], [yl(2)-0.1 yl(2)-0.1 yl(2)-0.05 yl(2)-0.05],'y')
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';

        alpha(0.3);
        str = sprintf('*');
        text(ntim(i)-3,yl(2)-0.03,str,'FontSize',20);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

 %% (3)load graphs, edge = midVOT
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
% remove the negative values after bias correction
for ipat = 1:np
    for iwin = 1:nwin
        mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1};
        mi456(mi456<0)=0;
        miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1} = mi456;
    end
end

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
mid_VOT_energy = zeros(np,nwin);

for iwin = 1:41
   % figure;hold on;
    for ipat = 1:np
        
        mat = squeeze(diff_midVOT(ipat,:,:,iwin));
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
title('Energy of matrices: (MidVOT- BS) >0');
%title('Energy of matrices: (LowVOT- BS) <0');

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


%% combined midVOT figure
%let's plot the energy of baseline normalized graphs
figure('Position',[100 100 750 500]);
for i = 1:np-1
    if i~=1
        a=plot(ntim,mid_VOT_energy(i,:),'b-*','LineWidth',4);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    else
        plot(ntim,mid_VOT_energy(i,:),'b-*','LineWidth',4);
    end
    hold on;
end
plot(ntim,mid_VOT_energy(np,:),'r-*','LineWidth',4);
set(gca,'FontSize',16);
xlim([-25,425]);
ylim = yl;
%legend('Healthy subject','Aphasic subject','Location','bestoutside');
ylabel('Frobenius norm','FontSize',20);
xlabel('Time (ms)','FontSize',20);

%title('Energy of matrices: (MidVOT- BS) >0');      

set(gcf,'color','w');
box off;
yl = ylim;

for i = 1:length(ntim)
    patstats = stats_mid_gt_bs(22,:);
    if(patstats(i)==1)
        a=patch([ntim(i)-0.5*10 ntim(i)+0.5*10 ntim(i)+0.5*10 ntim(i)-0.5*10], [yl(2)-0.1 yl(2)-0.1 yl(2)-0.05 yl(2)-0.05],'y')
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';

        alpha(0.3);
        str = sprintf('*');
        text(ntim(i)-3,yl(2)-0.03,str,'FontSize',20);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end


%% let's try plotting this using the same clustering technique

%%  lowVOT case

healthy_matricesLOW = squeeze(sum(diff_lowVOT(1:np-1,:,:,:)));
healthy_matricesHIGH = squeeze(sum(diff_highVOT(1:np-1,:,:,:)));
healthy_matricesMID = squeeze(sum(diff_midVOT(1:np-1,:,:,:)));

aphasic_matrixLOW = squeeze(diff_lowVOT(end,:,:,:));
aphasic_matrixHIGH = squeeze(diff_highVOT(end,:,:,:));
aphasic_matrixMID = squeeze(diff_midVOT(end,:,:,:));

new_edge_matrices = aphasic_matrixLOW - healthy_matricesLOW;
new_edge_matrices(new_edge_matrices<0)=0;

M.matrices = new_edge_matrices;

mat = (squeeze(M.matrices(:,:,1)));
n = numel(mat(:,1));
gamma = 0.9;
omega = 1;
T=length(squeeze(M.matrices(1,1,:)));
B=spalloc(n*T,n*T,n*n*T+2*n*T);
twomu=0;
for s=1:T
    k=sum(squeeze(M.matrices(:,:,s)));
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:n]+(s-1)*n;
    B(indx,indx)=squeeze(M.matrices(:,:,s))-gamma*k'*k/twom;
end
twomu=twomu+2*omega*n*(T-1);
B = B + omega*spdiags(ones(n*T,2),[-n,n],n*T,n*T);
[S,Q] = genlouvain(B);
Q = Q/twomu;
partition_multi = reshape(S,n,T);

% What does the partition of nodes into communities look like?
% example of a partition with time
figure; imagesc('XData',[ntim],'YData',1:26,'CData',partition_multi);
ylim([1,26]);xlim([0,400]);
box off;
set(gcf,'color','w');
xlabel('Time (ms)','FontSize',16);
ylabel('Channel number','FontSize',16);
%colorbar('Ticks',[1,2,3]);
ax=gca;
ax.FontSize = 16;


ch = partition_multi(:,1);

%% let's repeat this partitioning multiple times, and save them

nperm = 500;
partition_indices = zeros(length(chlist),nperm);
for iperm = 1:nperm
    M.matrices = new_edge_matrices;
    mat = (squeeze(M.matrices(:,:,1)));
    n = numel(mat(:,1));
    gamma = 0.9;
    omega = 1;
    T=length(squeeze(M.matrices(1,1,:)));
    B=spalloc(n*T,n*T,n*n*T+2*n*T);
    twomu=0;
    for s=1:T
        k=sum(squeeze(M.matrices(:,:,s)));
        twom=sum(k);
        twomu=twomu+twom;
        indx=[1:n]+(s-1)*n;
        B(indx,indx)=squeeze(M.matrices(:,:,s))-gamma*k'*k/twom;
    end
    twomu=twomu+2*omega*n*(T-1);
    B = B + omega*spdiags(ones(n*T,2),[-n,n],n*T,n*T);
    [S,Q] = genlouvain(B);
    Q = Q/twomu;
    partition_multi = reshape(S,n,T); 
    partition_indices(:,iperm) = partition_multi(:,1);    
end

%% consensus clustering
figure('Position',[100 100 1500 500]); 
imagesc(partition_indices);
box off;
set(gcf,'color','w');
xlabel('Repetitions','FontSize',16);
ylabel('Channel number','FontSize',16);
colorbar('Ticks',[1,2,3,4]);
ax=gca;
ax.FontSize = 16;

D= agreement(partition_indices);
figure; imagesc(D)
ciu = consensus_und(D,1,1000);
ch =ciu;
list_ch = 1:length(chlist);
[~,ind] = sort(ch);

figure; imagesc(D(list_ch(ind),list_ch(ind)))
box off;
set(gcf,'color','w');
xlabel('Channels','FontSize',16);
ylabel('Channels','FontSize',16);
axis off;
axis square;
colorbar;
ax=gca;
ax.FontSize = 16;
%% lets plot the energies for these individual partitions

ch=ciu;
nbrainRegions = 4;

nearestNeighbors = cell(1,nbrainRegions);

nearestNeighbors{1} = find(ch==1);
nearestNeighbors{2} = find(ch==2);
nearestNeighbors{3} = find(ch==3);
nearestNeighbors{4} = find(ch==4);

sig_windows_aphasia = zeros(nbrainRegions,nwin);
for inode = 1:nbrainRegions
           
    low_VOT_energy = zeros(np,nwin);

                for iwin = 1:41
                   % figure;hold on;
                    for ipat = 1:np

                        mat = squeeze(diff_lowVOT(ipat,(nearestNeighbors{inode}),(nearestNeighbors{inode}),iwin));
                     %   bar(ipat,norm(mat,'fro'));
                        low_VOT_energy(ipat,iwin) = norm(mat,'fro');
                    end
                end

                %let's plot the energy of baseline normalized graphs
                figure;
                for i = 1:np-1
                    if i~=1
                        a=plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
                        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    else
                        plot(ntim,low_VOT_energy(i,:),'b-*','LineWidth',4);
                    end
                    hold on;
                end
                plot(ntim,low_VOT_energy(np,:),'r-*','LineWidth',4);


                legend('Healthy subject','Aphasic subject');
                ylabel('Frobenius norm');
                xlabel('Time');
                title('Energy of matrices: (LowVOT- BS) >0');

                axis square;
                set(gcf,'color','w');
                box off;

    X1 = low_VOT_energy;%

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

    
    
end

     %% plot this in the graphs (cluster 1)

load('ChanNamesandLocs.mat');
chanNamesLocsShort = chanNamesLocs(chlist,:);
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);
inode = 4;
chl = nearestNeighbors{inode};

empty_graph = zeros(length(chlist),length(chlist));
graph1 = zeros(length(chl),length(chl));

% load the mat file
M.matrices = new_edge_matrices;

mat = (mean(M.matrices(:,:,:),3)); % mean edge metric

gf=figure;
empty_graph(chl,chl) = mat(chl,chl);

G = graph(empty_graph,chnames);
G1 = graph(graph1,chnames(chl));

LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
pe=plot(G,'XData',cell2mat(chanNamesLocsShort(:,2)),'YData',cell2mat(chanNamesLocsShort(:,3)),'NodeFontSize',14,'LineWidth',LWidths); hold on;
p=plot(G1,'XData',cell2mat(chanNamesLocsShort(chl,2)),'YData',cell2mat(chanNamesLocsShort(chl,3)),'NodeFontSize',14);


numlinks = degree(G);
numlinks(numlinks>0) = round(sum(mat(chl,chl)).*20);
highlight(pe,numlinks>0,'NodeColor','r')
pe.MarkerSize = [4+numlinks];

p.NodeColor='r';
az=-90;
ele = 90;
view([az ele])
box off;
set(gcf,'color','w');
annotation(gf,'ellipse',...
    [0.176,0.166666666666667,0.722214285714286,0.721428571428572]);

annotation(gf,'line',[0.473214285714286 0.516071428571428],...
    [0.880952380952381 0.933333333333334]);

annotation(gf,'line',[0.558928571428571 0.516071428571428],...
    [0.885714285714286 0.933333333333334]);
axis off;

