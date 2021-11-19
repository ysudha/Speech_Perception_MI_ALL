%code to test adhd within subjects example, from online stats book.

adhd = [57    48    64    62
        27    42    48    49
        32    39    34    30
        31    23    25    34
        34    36    42    38
        38    34    40    36
        71    61    59    77
        33    36    42    51
        34    41    36    45
        53    57    67    42
        36    32    29    43
        42    48    50    57
        26    48    44    36
        52    39    57    58
        36    36    36    35
        55    48    57    60
        36    35    30    33
        42    39    36    49
        36    32    41    33
        54    48    65    59
        34    35    35    35
        29    28    40    37
        33    27    43    45
        33    40    44    29];
    
    [p,tbl,stats] = anova2(adhd)
    
    % 24 subjects, 4 dosages, is there a difference in the dosage
    % condition?

% ids of subjects
adhdIds =  cellstr(num2str([1:24]'));

t = table(adhdIds,adhd(:,1),adhd(:,2),adhd(:,3),adhd(:,4),'VariableNames',{'IDs','d1','d2','d3','d4'});

Meas = dataset([1 2 3 4]','VarNames',{'Measurements'});

rm = fitrm(t,'d1-d4 ~ IDs-1','WithinDesign',Meas);

ranovatbl = ranova(rm)





% %% make it horizontal
% 
% adhdIds =  cellstr(strcat('ID',(num2str([10:33]'))));
% adhdIds = adhdIds';
% 
% dosageIds = {'d1';'d2';'d3';'d4'};
% 
% t = table(dosageIds,adhd(1,:)',adhd(2,:)',adhd(3,:)',adhd(4,:)',adhd(5,:)'...
%     ,adhd(6,:)',adhd(7,:)',adhd(8,:)',adhd(9,:)',adhd(10,:)',adhd(11,:)',...
%     adhd(12,:)',adhd(13,:)',adhd(14,:)',adhd(15,:)',adhd(16,:)',adhd(17,:)',...
%     adhd(18,:)',adhd(19,:)',adhd(20,:)',adhd(21,:)',adhd(22,:)',adhd(23,:)',...
%     adhd(24,:)','VariableNames',{'dosIDs',adhdIds{1:24}});
% 
%     
% Meas = dataset([1:24]','VarNames',{'Measurements'});
% 
% rm = fitrm(t,'ID10-ID33 ~ dosIDs-1','WithinDesign',Meas);
% 
% ranovatbl = ranova(rm)