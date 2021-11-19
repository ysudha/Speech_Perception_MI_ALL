% Example datasets that work

%% Crawford & Howell's modified t-test
%  Crawford JR, Garthwaite PH. Neuropsychology 2005 May;19(3):318-31)
Xm = 60;
Xs = 7;
Xc = 33;
Ym = 24;
Ys = 4.8;
Yc = 15;
n = 20;
r = 0.68;
alpha = 0.05;
[h p t df] = ttestch(60, 7, 33, 20, 0.05)
% should be [1 <0.001 -3.7642 19]
[h p t df] = ttestch(24, 4.8, 15, 20, 0.05)
% should be [1 0.042 -1.8298 19]

%% http://slideplayer.com/slide/778140/
Xm = 100;
Xs = 14.32;
Xc = 74;
Ym = 50;
Ys = 12.34;
Yc = 61;
n = 20;
r = 0.72;
alpha = 0.05;
[h p t df] = ttestch(100, 14.32, 74, 20, 0.05)
% should be [1 0.0462 -1.7719 19]
[h p t df] = ttestch(50, 12.34, 61, 20, 0.05)
% should be [1 0.1979 0.8699 19]
