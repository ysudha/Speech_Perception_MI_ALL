function [MI] = MutualInformation(x,y,K)
%entropy of x and y jointly.
if (size(x,1) < size(x,2))
    x = x.';
end
if (size(y,1) < size(y,2))
    y = y.';
end

dataMat = [x y];
NumPoints = size(y,1);
dim_xy = size(dataMat,2);
dim_x = size(x,2);
dim_y = size(y,2);

atria_xy = nn_prepare(dataMat,'maximum');
[indexMat, distance_xy] = nn_search(dataMat,atria_xy,1:NumPoints,K,0);

%H_XY = -psi(K) + psi(NumPoints) + (dim_xy)*mean(log(2*distance_xy(:,K)));

%entropy of x
atria_x = nn_prepare(x,'maximum');
ncnt_x = range_search(x,atria_x,1:NumPoints,distance_xy(1:NumPoints,K)-eps,0);
%H_X = -mean(psi(ncnt_x+1)) + psi(NumPoints) + (dim_x)*mean(log(2*distance_xy(:,K)));

%entropy of y
atria_y = nn_prepare(y,'maximum');
ncnt_y = range_search(y,atria_y,1:NumPoints,distance_xy(1:NumPoints,K)-eps,0);
%H_Y = -mean(psi(ncnt_y+1)) + psi(NumPoints) + (dim_y)*mean(log(2*distance_xy(:,K)));


%MI KSG 1
%MI = H_X+H_Y-H_XY
MI = psi(K)  - mean(psi(ncnt_x+1) + psi(ncnt_y+1)) + psi(NumPoints);
end