% This one is to plot the raw UNNORMALIZED forbenius norm of the adjacency
% matrices, including baseline.


%% plot raw MI matrices

 %% (1)  load graphs, edge = lowVOT
chlist = [1:25 28];
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = ([1 51 101:5:301]*2-202) + 50;
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

%Why? Leys look at the absolute value of normalized matrices themselves
alpha_val=0.5;
figure; hold on;
for ipat = 1:22
    %figure;
    vec_frob_norm = zeros(1,nwin);
    for iwin = 1:43 
        mat = squeeze(miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1});
        vec_frob_norm(1,iwin) = norm(mat,'fro');
        
    end
    if(ipat<22)
        plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[27,158,119,alpha_val*255]/255, 'MarkerSize',2);
    else
        plot(ntim,vec_frob_norm,'r-o','LineWidth',3,'MarkerSize',3);
    end
end
set(gca,'FontSize',16);
set(gcf,'color','w');
box off;

%repeat figure with legend
figure; hold on;
for ipat = 1:22
    %figure;
    vec_frob_norm = zeros(1,nwin);
    for iwin = 1:43 
        mat = squeeze(miVOT123.mi_matrix_VOT123_K4{ipat,1}{iwin,1});
        vec_frob_norm(1,iwin) = norm(mat,'fro');
        
    end
    if(ipat<22 & ipat >1)
        a=plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[27,158,119,alpha_val*255]/255, 'MarkerSize',2);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    elseif(ipat==1)
        plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[27,158,119,alpha_val*255]/255, 'MarkerSize',2);
    else
        plot(ntim,vec_frob_norm,'r-o','LineWidth',2,'MarkerSize',2);
    end
end
set(gca,'FontSize',16);
set(gcf,'color','w');
box off;
legend('Healthy subjects LowVOT','Aphasic subject LowVOT','Location','bestoutside');
legend box off;
%% (2)   load graphs, edge = midVOT
chlist = [1:25 28];
miVOT456 = load('mi_matrix_VOT456_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = ([1 51 101:5:301]*2-202) + 50;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT

for ipat = 1:np
    for iwin = 1:nwin
        mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1};
        mi456(mi456<0)=0;
        miVOT456.mi_matrix_VOT123_K4{ipat,1}{iwin,1} = mi456;
    end
end
%Why? Leys look at the absolute value of normalized matrices themselves
alpha_val=0.5;
figure; hold on;
for ipat = 1:22
    %figure;
    vec_frob_norm = zeros(1,nwin);
    for iwin = 1:43 
        mat = squeeze(miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1});
        vec_frob_norm(1,iwin) = norm(mat,'fro');
        
    end
    if(ipat<22)
        plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[217,95,2,alpha_val*255]/255, 'MarkerSize',2);
    else
        plot(ntim,vec_frob_norm,'r-o','LineWidth',3,'MarkerSize',3);
    end
end
set(gca,'FontSize',16);
set(gcf,'color','w');
box off;

%repeat figure with legend
figure; hold on;
for ipat = 1:22
    %figure;
    vec_frob_norm = zeros(1,nwin);
    for iwin = 1:43 
        mat = squeeze(miVOT456.mi_matrix_VOT456_K4{ipat,1}{iwin,1});
        vec_frob_norm(1,iwin) = norm(mat,'fro');
        
    end
    if(ipat<22 & ipat >1)
        a=plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[217,95,2,alpha_val*255]/255, 'MarkerSize',2);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    elseif(ipat==1)
        plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[217,95,2,alpha_val*255]/255, 'MarkerSize',2);
    else
        plot(ntim,vec_frob_norm,'r-o','LineWidth',2,'MarkerSize',2);
    end
end
set(gca,'FontSize',16);
set(gcf,'color','w');
box off;
legend('Healthy subjects MidVOT','Aphasic subject MidVOT','Location','bestoutside');
legend box off;
%% (3)   load graphs, edge = highVOT
chlist = [1:25 28];
miVOT789 = load('mi_matrix_VOT789_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = length(chlist);
ntim = ([1 51 101:5:301]*2-202) + 50;
% creating the adjacency matrices, with edge = VOT123, only
% positive increases diff_mid_lowVOT
diff_highVOT = zeros(np,nch,nch,nwin-2);
alpha_val=0.3;
for ipat = 1:np
    for iwin = 1:nwin
        mi789 = miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1};
        mi789(mi789<0)=0;
        miVOT789.mi_matrix_VOT123_K4{ipat,1}{iwin,1} = mi789;
    end
end
figure; hold on;
for ipat = 1:22
    %figure;
    vec_frob_norm = zeros(1,nwin);
    for iwin = 1:43 
        mat = squeeze(miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1});
        vec_frob_norm(1,iwin) = norm(mat,'fro');
        
    end
    if(ipat<22)
        plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[117,112,179,alpha_val*255]/255, 'MarkerSize',2);
    else
        plot(ntim,vec_frob_norm,'r-o','LineWidth',2,'MarkerSize',2);
    end
end
set(gca,'FontSize',16);
set(gcf,'color','w');
box off;

%repeat figure with legend
figure; hold on;
for ipat = 1:22
    %figure;
    vec_frob_norm = zeros(1,nwin);
    for iwin = 1:43 
        mat = squeeze(miVOT789.mi_matrix_VOT789_K4{ipat,1}{iwin,1});
        vec_frob_norm(1,iwin) = norm(mat,'fro');
        
    end
    if(ipat<22 & ipat >1)
        a=plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[117,112,179,alpha_val*255]/255, 'MarkerSize',2);
        a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    elseif(ipat==1)
        plot(ntim,vec_frob_norm,'o-','LineWidth',2,'Color',[117,112,179,alpha_val*255]/255, 'MarkerSize',2);
    else
        plot(ntim,vec_frob_norm,'r-o','LineWidth',2,'MarkerSize',2);
    end
end
set(gca,'FontSize',16);
set(gcf,'color','w');
box off;
legend('Healthy subjects HighVOT','Aphasic subject HighVOT','Location','bestoutside');
legend box off;