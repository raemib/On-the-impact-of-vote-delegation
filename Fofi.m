%% TEST --- Conventional Voting VS Free Vote Delegation
clear

P=load('Pascal.mat');
P=P.P;

%P(i,j) = n cho 

cap = 2;

mini = 3;
maxi = 3;

minm = 1;
maxm = 1;

pmin = 0.7850;
pmax = 0.7950;
pstep = 0.0000001;

%p=0.60; % fix p for Mesh Plot

%Result = zeros(maxi,7);
j = 1;

for p = pmin:pstep:pmax
    D1=0;
    for t = mini:maxi
        for m = minm:maxm
            %F = Fofi(maxi,k);
            CV = ConVoting(p,t,P);
            FVD = FreeVoteDeleg(p,t,m,P);
            %CVD = CappedVoteDeleg(p,t,m,cap,P);
            Delta = CV - FVD; % change for FVD or CVD
            if Delta > D1
                MAXD = Delta;
                MAXCV = CV;
                MAXFVD = FVD; % change for FVD or CVD
                MAXt = t;
                MAXi = maxi;
                MAXm = m;
                MAXp = p;
                D1 = Delta;
            end
            % for Mesh Plot
            %M(m,t) = m;
            %K(m,t) = t;
            %D(t,m) = Delta;
        end
    end
    Result(j,1) = MAXD;
    Result(j,2) = MAXCV;
    Result(j,3) = MAXFVD;
    Result(j,4) = MAXt;
    Result(j,5) = MAXi;
    Result(j,6) = MAXm;
    Result(j,7) = MAXp;
    j = j+1;
end

%M = M.';

T = table(Result(:,1),Result(:,2),Result(:,3),Result(:,4),Result(:,5),Result(:,6),Result(:,7));
T.Properties.VariableNames = {'Delta','CV','FVD','t','max i','m','p'};


clearvars -except T K M D

plot(table2array(T(:,7)),table2array(T(:,1)))

%% Conventional Voting, Free Delegation, Capped Delegation
clear

P=load('Pascal.mat');
P=P.P;

cap = 2;

maxi = 20;

minm = 1;
maxm = 20;

pmin = 0.7;
pmax = 0.7;
pstep = 0.05;

%p=0.60; % fix p for Mesh Plot

%Result = zeros(maxi,7);
j = 1;

D1 = 0;
for p = pmin:pstep:pmax
    for t = 1:maxi
        for m = minm:maxm
            %F = Fofi(maxi,k);
            CV = ConVoting(p,t,P);
            FVD = FreeVoteDeleg(p,t,m,P);
            CVD = CappedVoteDeleg(p,t,m,cap,P);
            DeltaCVFVD = CV - FVD;
            DeltaCVCVD = CV - CVD;
            %if DeltaFVD > D1
                Result(j,1) = CV;
                Result(j,2) = FVD;
                Result(j,3) = CVD;
                Result(j,4) = DeltaCVFVD;
                Result(j,5) = DeltaCVCVD;
                Result(j,6) = t;
                Result(j,7) = maxi;
                Result(j,8) = m;
                Result(j,9) = p;
                j = j+1;
                %D1 = DeltaFVD;
            %end
            % for Mesh Plot
            %M(m,k) = m;
            %K(k,m) = k;
            %D(k,m) = Delta;
        end
    end
end

%M = M.';

T = table(Result(:,1),Result(:,2),Result(:,3),Result(:,4),Result(:,5),Result(:,6),Result(:,7),Result(:,8),Result(:,9));
T.Properties.VariableNames = {'CV','FVD','CVD','CV-FVD','CV-CVD','k','max i','m','p'};

clearvars -except T K M D

%% Conventional Voting VS Free Vote Delegation
clear

P=load('Pascal.mat');
P=P.P;

%P(i,j) = n cho 

mini = 1;
maxi = 50;

minm = 1;
maxm = 50;

pmin = 0.6;
pmax = 0.6;
pstep = 0.00001;

%p=0.60; % fix p for Mesh Plot

%Result = zeros(maxi,7);
j = 1;

D1 = 0;
for p = pmin:pstep:pmax
    for t = mini:maxi
        for m = minm:maxm
            %F = Fofi(maxi,k);
            CV = ConVoting(p,t,P);
            FVD = FreeVoteDeleg(p,t,m,P);
            Delta = CV - FVD;
            %if Delta > D1
                Result(j,1) = Delta;
                Result(j,2) = CV;
                Result(j,3) = FVD;
                Result(j,4) = t;
                Result(j,5) = maxi;
                Result(j,6) = m;
                Result(j,7) = p;
                j = j+1;
                %D1 = Delta;
            %end
            % for Mesh Plot
            M(m,t) = m;
            K(m,t) = t;
            D(t,m) = Delta;
        end
    end
end


T = table(Result(:,1),Result(:,2),Result(:,3),Result(:,4),Result(:,5),Result(:,6),Result(:,7));
T.Properties.VariableNames = {'Delta','CV','FVD','t','max i','m','p'};


clearvars -except T K M D p

%% 2D Plot
MP=M.';
for i = 1:7
    Diff(i,1:30) = D(i,1:30);
    Mdiff(i,1:30) = MP(i,1:30);
    Diff2(i,1:30) = D(1:30,i);
    Tdiff(i,1:30) = K(i,1:30);
end

figure(1)
for i = 1:7
    hold on
    plot(Mdiff(i,:),Diff(i,:),'LineWidth',1.5)
end
legend('m=1','m=2','m=3','m=4','m=5','m=6','m=7')

xlabel('t')
ylabel('P(p) - P(p,m)')
set(gca,'fontsize', 18)
title(['p = ',num2str(p)], 'Fontsize', 24)

hold off


figure(2)
for i = 1:7
    hold on
    plot(Tdiff(i,:),Diff(i,:),'LineWidth',1.5)
end
legend('m=1','m=2','m=3','m=4','m=5','m=6','m=7')

xlabel('t')
ylabel('P(p) - P(p,m)')
set(gca,'fontsize', 18)
title(['p = ',num2str(p)], 'Fontsize', 24)

%% Conventional Voting VS Capped Vote Delegation
%clear

P=load('Pascal.mat');
P=P.P;

cap = 2;

mini = 1;
maxi = 10;

minm = 1;
maxm = 10;

pmin = 0.70;
pmax = 0.70;
pstep = 0.05;

p=0.60; % fix p for Mesh Plot

%Result = zeros(maxi,7);
j = 1;

D1 = 0;
for p = pmin:pstep:pmax
    for t = mini:maxi
        for m = minm:maxm
           % F = Fofi(maxi,k);
            CV = ConVoting(p,t,P);
            CVD = CappedVoteDeleg(p,t,m,cap,P);
            Delta = CV - CVD;
            if Delta > D1
                Result(j,1) = Delta;
                Result(j,2) = CV;
                Result(j,3) = CVD;
                Result(j,4) = t;
                %Result(j,5) = maxi;
                Result(j,5) = m;
                Result(j,6) = p;
                j = j+1;
                D1 = Delta;
            end
            % for Mesh Plot
            M(m,t) = m;
            K(m,t) = t;
            D(t,m) = Delta;
        end
    end
end


T = table(Result(:,1),Result(:,2),Result(:,3),Result(:,4),Result(:,5),Result(:,6));
T.Properties.VariableNames = {'Delta','CV','CVD','t','m','p'};

clearvars -except T K M D SAVED

%% Mesh Plot
figure(1)
mesh(Mnew,Knew,Dnew)
xlabel('t')
ylabel('m')
zlh = zlabel('P(p) - Pc(p,m)');
%zlh.Position(1) = -16.6855;
%zlh.Position(2) = 109.8340;
%zlh.Position(3) = 0.0110;
xticks([0  20  40  60  80  100])
yticks([0  20  40  60  80  100])
zticks([0  0.025  0.050 ])
title(['p = ',num2str(p)], 'Fontsize', 24)
zlim([0 0.05])
pbaspect([1 1 0.5])
view(-45,20)
set(gca,'fontsize', 18)


%% F(i) = 0, F(k) = 1
function [F] = Fofi(maxi,k)
F = zeros(maxi,1);
F(k)=1;
end

%% CV
function [CV] = ConVoting(p,t,P)
CV = 0;
i = t;
value = 0;
   for k=0:i
       value = value + P(i-k+1,k+1) * p^k * (1-p)^(i-k) * gfunction(k,i-k);
       % P(i-k+1,k+1) = i choose k
   end
CV = value;
end

%% FVD
function [FVD] = FreeVoteDeleg(p,t,m,P)
FVD = 0;
i = t;
    value = 0;
    for k = 0:i
        value = value + P(i-k+1,k+1) * p^k * (1-p)^(i-k) * Gfunction(k,i-k,m,P);
        % P(i-k+1,k+1) = i choose k
    end
    FVD = value;
end

%% CVD
function [CVD] = CappedVoteDeleg(p,t,m,cap,P)
CVD = 0;
i=t;
value = 0;
for k = 0:i
    value = value + P(i-k+1,k+1) * p^k * (1-p)^(i-k) * Gcapfunction(k,i-k,m,cap,P);
    % P(i-k+1,k+1) = i choose k
CVD = value;
end
end

%% G function Capped
function [G] = Gcapfunction(k,l,m,cap,P)
G = 0;
if k == l
    G = 0.5;
elseif k == 0
    G = 0;
elseif l == 0
    G = 1;
else
    for h = 0:m
        G = G + P(m-h+1,h+1) * (k/(k+l))^h * (l/(k+l))^(m-h) * gfunction(min(k+h,cap*k),min(l+m-h,cap*l));
        % P(m-h+1,h+1) = m choose h
    end
end
end

%% G function
function [G] = Gfunction(k,l,m,P)
G = 0;
if k == l
    G = 0.5;
elseif k == 0
    G = 0;
elseif l == 0
    G = 1;
else
    for h = 0:m
        G = G + P(m-h+1,h+1) * (k/(k+l))^h * (l/(k+l))^(m-h) * gfunction(k+h,l+m-h);
        % P(m-h+1,h+1) = m choose h
    end
end
end

%% g function
function [g] = gfunction(k,l)
if k>l
    g = 1;
elseif k==l
    g = 0.5;
else
    g = 0;
end
end

%% Create the Pascal Matrix P
function[P] = CalcPascal()
% Output: P(i-j+1,j+1) = i choose j
% Define size of the Matrix and initialize it
Size = 10000;
P = zeros(Size, Size);
% Create the first two columns
P(:,1) = ones(Size,1);
a = [1:1:Size]';
P(:,2) = a;
P(1,:) = ones(1,Size);

% Create the rest of the matrix
for i = 2:Size
    for j = 3:Size
        P(i,j) = P(i-1,j) + P(i,j-1);
    end
end
end


