function [tstats_sorted,thresh] = func_cluster_analysis(nperm,nwin,np,low_VOT,mid_VOT,high_VOT,max_t_cluster)
%nperm = 500;
max_tstats = zeros(1,nperm);
for iperm = 1:nperm
    pval_mid_gt_low_and_high = zeros(1,nwin);
    contrasts = zeros(np-1,nwin);
    t_statistics = zeros(1,nwin);
    for i = 1:nwin
        anovaMatrix = [low_VOT(:,i) mid_VOT(:,i)  high_VOT(:,i)];
        bb = normalize_multiple_patients(anovaMatrix);
        for j = 1:np-1
            bb(j,:) = bb(j,randperm(3));
        end
        %bb = anovaMatrix;
        col = bb(:,2) - 0.5*(bb(:,1)+bb(:,3));
        contrasts(:,i) = col;
        [h,p,ci,stats] = ttest(col,0,'Tail','right');
        pval_mid_gt_low_and_high(i)=p;
        t_statistics(i) = stats.tstat;
    end

tmask = t_statistics>1.725; % for df=21-1, critical t stat for 95% confidence is 1.725
cluster_starts = strfind((t_statistics>1.725),[0 1])+1; %+1 as we dont want to include 0 in the string 011
if(tmask(1)==1)
    cluster_starts = [1 cluster_starts];
end
t_stats_clusters = zeros(1,length(cluster_starts));
for j = 1:length(cluster_starts)
    
    within_cluster_incementer = cluster_starts(j);
    while(tmask(within_cluster_incementer)==1 && within_cluster_incementer<41)
        t_stats_clusters(j) = t_stats_clusters(j)+ t_statistics(within_cluster_incementer);
        within_cluster_incementer = within_cluster_incementer+1;
    end
end
  
if(length(t_stats_clusters)>0)
    max_tstats(iperm) = max(t_stats_clusters);
else
    max_tstats(iperm) = 0;
end

end

tstats_sorted = sort(max_tstats);
thresh = tstats_sorted(round(nperm*0.95))

%pval = sum(max_tstats>max_t_cluster)./nperm;
