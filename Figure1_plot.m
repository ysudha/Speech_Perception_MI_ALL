% Figure with empty graph, and raw MI matrices (presently fig2)

%% plot the empty graph for figure 
chlist = [1:25 28];

load('ChanNamesandLocs.mat');
chnamesOrig = chanNamesLocs(:,1);
chnames = chnamesOrig(chlist);

empty_graph = zeros(length(chlist),length(chlist));


G = graph(empty_graph,chnames);
gf=figure;
p=plot(G,'XData',cell2mat(chanNamesLocs(chlist,2)),'YData',cell2mat(chanNamesLocs(chlist,3)),'NodeFontSize',14); hold on;
p.NodeColor='r';
hold on;
az=-90;
ele = 90;
view([az ele])
box off;
axis off;
set(gcf,'color','w');
annotation(gf,'ellipse',...
    [0.176,0.166666666666667,0.722214285714286,0.721428571428572]);

annotation(gf,'line',[0.473214285714286 0.516071428571428],...
    [0.880952380952381 0.933333333333334]);

annotation(gf,'line',[0.558928571428571 0.516071428571428],...
    [0.885714285714286 0.933333333333334]);



%% plot raw MI matrices

 %% (1)  load graphs, edge = lowVOT
chlist = [1:25 28];
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = [1 51 101:5:301];
ntim = [101:5:301]*2-202;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases 
diff_lowVOT = zeros(np,nch,nch,nwin-2);

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


%Why? Leys look at the absolute value of normalized matrices themselves
for ipat = [22]
    %figure;
    index=1;
    for iwin = [1 2 39 40 41] 
        mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
        %subplot(1,3,index);
        figure;
        imagesc(mat);
        colormap jet; 
        %colorbar('fontweight','bold', 'fontSize',18);
        caxis([0 0.15]);
        axis square;
        index=index+1;
        set(gcf,'color','w');
        axis off;
    end
end

%% trial windows and baseline RAW MI - just 3 windows

% remove the negative values after permuation correction
for ipat = 1:np
    for iwin = 1:nwin
        mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1};
        mi123(mi123<0)=0;
        miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1} = mi123;
    end
end

for ipat = 1
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    for iwin = [1 3 43]
       mat2 = squeeze(mi123{iwin,1}); 
       mat = squeeze(mat2);
       mat = mat(chlist,chlist);
       
        figure;
        imagesc(mat);
        colormap jet;% colorbar;
        colorbar('fontweight','bold', 'fontSize',18);
        caxis([0 1]);
        axis square;
        set(gcf,'color','w');
        axis off;
    end

end

%% trial windows &  baseline - many windows

for ipat = [1]
        mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
for iwin = 1:43
       mat2 = squeeze(mi123{iwin,1}); 
       mat = squeeze(mat2);
       mat = mat(chlist,chlist);
       
        figure;
        imagesc(mat);
        colormap jet;% colorbar;
        %colorbar('fontweight','bold', 'fontSize',18);
        caxis([0 1]);
        axis square;
        set(gcf,'color','w');
        axis off;
    end
end


%% baseline normalized trial windows - many windows

for ipat = [1]
            mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows

for iwin = 3:43
       mat2 = squeeze(mi123{iwin,1} - mi123{1,1}); 
       mat = squeeze(mat2);
       mat = mat(chlist,chlist);
       
        figure;
        imagesc(mat);
        colormap jet; colorbar;
        colorbar('fontweight','bold', 'fontSize',18);
        caxis([-0.1 0.1]);
        axis square;
        set(gcf,'color','w');
        axis off;
    end
end

%% absolute value of baseline normalized graphs 

for ipat = [1]
    for iwin = 1:nwin 
        mat = squeeze(diff_lowVOT(ipat,:,:,iwin));
        figure;
        imagesc(mat);
        colormap jet; 
        %colorbar('fontweight','bold', 'fontSize',18);
        caxis([0 0.15]);
        axis square;
        set(gcf,'color','w');
        axis off;
    end
end

%% (2)   load graphs, edge = midVOT
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

for ipat = 1:np
    for iwin = 1:nwin
        mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1};
        mi456(mi456<0)=0;
        miVOT456.mi_matrix_VOT123_K4{ipat,1}{iwin,1} = mi456;
    end
end

for ipat = 1:np
    mi456 = miVOT456.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    for iwin = 3:nwin
       mat2 = squeeze(mi456{iwin,1} - mi456{1,1}); 
       mat = squeeze(mat2);
       diff_midVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

diff_midVOT(diff_midVOT<0)=0; %< is the right one
diff_midVOT = abs(diff_midVOT);


%Why? Leys look at the absolute value of normalized matrices themselves
for ipat = [22]
    %figure;
    index=1;
    for iwin = [1 2 39 40 41] 
        mat = squeeze(diff_midVOT(ipat,:,:,iwin));
        %subplot(1,3,index);
        figure;
        imagesc(mat);
        colormap jet; 
        %colorbar('fontweight','bold', 'fontSize',18);
        caxis([0 0.15]);
        axis square;
        index=index+1;
        set(gcf,'color','w');
        axis off;
    end
end


%% (3)   load graphs, edge = highVOT
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

for ipat = 1:np
    for iwin = 1:nwin
        mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1};
        mi789(mi789<0)=0;
        miVOT789.mi_matrix_VOT123_K4{ipat,1}{iwin,1} = mi789;
    end
end

for ipat = 1:np
    mi789 = miVOT789.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    for iwin = 3:nwin
       mat2 = squeeze(mi789{iwin,1} - mi789{1,1}); 
       mat = squeeze(mat2);
       diff_highVOT(ipat,:,:,iwin-2) = mat(chlist,chlist);
    end

end

nwin = nwin-2;

diff_highVOT(diff_highVOT<0)=0; %< is the right one
diff_highVOT = abs(diff_highVOT);


%Why? Leys look at the absolute value of normalized matrices themselves
for ipat = [22]
    %figure;
    index=1;
    for iwin = [1 2 39 40 41] 
        mat = squeeze(diff_highVOT(ipat,:,:,iwin));
        %subplot(1,3,index);
        figure;
        imagesc(mat);
        colormap jet; 
        %colorbar('fontweight','bold', 'fontSize',18);
        caxis([0 0.15]);
        axis square;
        index=index+1;
        set(gcf,'color','w');
        axis off;
    end
end
