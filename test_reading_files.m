% load the data of 21 healthy undergrads
%load('21SubjCell_ElecMatrixSeparatedfromTrialInfo.mat');

% data is in subjCell 21 x 2500 x 2
% lets see what is in column 2,3 4 etc of the trial id

test=[];
for i = 1:1900
    a = subjCell{1,i,2};
    test = [test a(5)];
end

