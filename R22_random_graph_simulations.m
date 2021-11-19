
%% single graph
n = 150;
p = 0.01;

rand('seed',100); % reseed so you get a similar picture
G = rand(n,n) < p;
G = triu(G,1);
G = G + G';
Gr = graph(G);

thresh = 1/(n-1) %thresh<p => big component exists
figure;
plot(Gr,'Layout','force')

N = numedges(Gr)

[U1,S1,V1] = eig(G);
energy_G = sum(abs(diag(S1)))

%% loop over inreasing p values

n = 150;
p = 0.00001:0.0001:0.01;
N = zeros(size(p));
energy_ER = zeros(size(p));
ind=1;
for p = 0.00001:0.0001:0.01

    %rand('seed',100); % reseed so you get a similar picture
    G = rand(n,n) < p;
    G = triu(G,1);
    G = G + G';
    Gr = graph(G);

    thresh = 1/(n-1); %thresh<p => big component exists
   % figure;
   % plot(Gr,'Layout','force')

    N(ind) = numedges(Gr);

    [U1,S1,V1] = eig(G);
    energy_ER(ind) = sum(abs(diag(S1)));
    ind = ind+1;
end

figure;
plot(0.00001:0.0001:0.01,energy_ER)
figure;
plot(N,energy_ER)
%% one fixed graph with a large componenet, then randomize it

n = 150;
p = 0.01;

%rand('seed',100); % reseed so you get a similar picture
G = rand(n,n) < p;
G = triu(G,1);
G = G + G';
Gr = graph(G);

thresh = 1/(n-1) %thresh<p => big component exists
figure;
plot(Gr,'Layout','force')

N = numedges(Gr)

[U1,S1,V1] = eig(G);
energy_G = sum(abs(diag(S1)))

for i = 1:1
    rp = randperm(n);
    D = G(rp,rp);
    figure;
    plot(graph(D),'Layout','force')
end