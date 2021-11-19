% Example datasets that work

%% Revised Standardized Difference Test
% Crawford JR, Garthwaite PH. Neuropsychology 2005 May;19(3):318-31)
% Two tailed test
Xm = 60;
Xs = 7;
Xc = 33;
Ym = 24;
Ys = 4.8;
Yc = 15;
n = 20;
r = 0.68;
alpha = 0.05;
[h p t df] = rsdt(60, 7, 33, 24, 4.8, 15, 20, 0.68, 0.05)
% should be [1 0.033 2.2998 19]

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
[h p t df] = rsdt(100, 14.32, 74, 50, 12.34, 61, 20, 0.72, 0.05)
% should be [1 0.0035 3.3356 19]