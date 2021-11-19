 %% function plot this in the graphs binary
 % optionVal = 1, plot only the first 6 nodes of high degrees
 %option Value not equal to 1, plot the graph as well

function[]= plot_on_brain(new_edge_matrices,small_subset_time,optionVal,nnodes)

if(optionVal==1)
    chlist = [1:25 28];
    % load the mat file
    M.matrices = new_edge_matrices(:,:,small_subset_time);

    mat = (mean(M.matrices(:,:,:),3)); % mean edge metric
    mat(mat>0)=1;

    deg_matrix = sum(mat);
    [val,ind] = sort(deg_matrix,'descend')

    load('ChanNamesandLocs.mat');
    chanNamesLocsShort = chanNamesLocs(chlist,:);
    chnamesOrig = chanNamesLocs(:,1);
    chnames = chnamesOrig(chlist);
    inode = ind(1:nnodes);
    chl = inode;

    empty_graph = zeros(length(chlist),length(chlist));
    graph1 = zeros(length(chl),length(chl));


    gf=figure;

    G = graph(empty_graph,chnames);
    G1 = graph(graph1,chnames(chl));

    pe=plot(G,'XData',cell2mat(chanNamesLocsShort(:,2)),'YData',cell2mat(chanNamesLocsShort(:,3)),'NodeFontSize',14); hold on;
    p=plot(G1,'XData',cell2mat(chanNamesLocsShort(chl,2)),'YData',cell2mat(chanNamesLocsShort(chl,3)),'NodeFontSize',14);

    numlinks = 2*nnodes:-2:1;
    %highlight(p,numlinks,'NodeColor','r')
    p.MarkerSize = [4+numlinks];

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
else
    chlist = [1:25 28];
    % load the mat file
    M.matrices = new_edge_matrices(:,:,small_subset_time);

    mat = (mean(M.matrices(:,:,:),3)); % mean edge metric
    mat(mat>0)=1;

    deg_matrix = sum(mat);
    [val,ind] = sort(deg_matrix,'descend')

    load('ChanNamesandLocs.mat');
    chanNamesLocsShort = chanNamesLocs(chlist,:);
    chnamesOrig = chanNamesLocs(:,1);
    chnames = chnamesOrig(chlist);
    inode = ind(1:nnodes);
    chl = inode;

    empty_graph = mat;%zeros(length(chlist),length(chlist));
    graph1 = mat(chl,chl);


    gf=figure;

    G = graph(empty_graph,chnames);
    G1 = graph(graph1,chnames(chl));

    pe=plot(G,'XData',cell2mat(chanNamesLocsShort(:,2)),'YData',cell2mat(chanNamesLocsShort(:,3)),'NodeFontSize',14); hold on;
    p=plot(G1,'XData',cell2mat(chanNamesLocsShort(chl,2)),'YData',cell2mat(chanNamesLocsShort(chl,3)),'NodeFontSize',14);

    numlinks = 2*nnodes:-2:1;
    %highlight(p,numlinks,'NodeColor','r')
    p.MarkerSize = [4+numlinks];

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
end
    
end