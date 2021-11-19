   conditions = cell(21*3,1);
    for k = 1:21
        conditions{k} = 'lowVOT';
    end
    for k = 22:22+20
        conditions{k} = 'midVOT';
    end
    for k = 43:43+20
        conditions{k} = 'highVOT';
    end
     t = table(conditions,[bb(:,1);bb(:,2);bb(:,3)],'VariableNames',{'conditions','p'});
    Meas = table([1]','VariableNames',{'Measurements'})
rm = fitrm(t,'p~conditions','WithinDesign',Meas)    
ranovatbl = ranova(rm,'WithinDesign',Meas)
[p,tbl, stats]=anova1(anovaMatrix)

