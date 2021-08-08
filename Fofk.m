%% Conventional Voting VS Free Vote Delegation

clear

maxi = 10;

minm = 1;
maxm = 20;

pmin = 0.55;
pmax = 0.95;
pstep = 0.05;

p=0.60; % fix p for Mesh Plot

%Result = zeros(maxi,7);
j = 1;

D1 = 0;
for p = pmin:pstep:pmax
    for k = 1:maxi
        for m = minm:maxm
            F = Fofi(maxi,k);
            CV = ConVoting(p,maxi,F);
            FVD = FreeVoteDeleg(p,maxi,F,m);
            Delta = CV - FVD;
            if Delta > D1
                Result(j,1) = Delta;
                Result(j,2) = CV;
                Result(j,3) = FVD;
                Result(j,4) = k;
                Result(j,5) = maxi;
                Result(j,6) = m;
                Result(j,7) = p;
                j = j+1;
                D1 = Delta;
            end
            % for Mesh Plot
            %M(m,k) = m;
            %K(k,m) = k;
            %D(k,m) = Delta;
        end
    end
end

%M = M.';

T = table(Result(:,1),Result(:,2),Result(:,3),Result(:,4),Result(:,5),Result(:,6),Result(:,7));
T.Properties.VariableNames = {'Delta','CV','FVD','k','max i','m','p'};

clearvars -except T K M D

%% Conventional Voting VS Capped Vote Delegation
clear

cap = 2;

maxi = 10;

minm = 1;
maxm = 20;

pmin = 0.55;
pmax = 0.95;
pstep = 0.05;

p=0.60; % fix p for Mesh Plot

%Result = zeros(maxi,7);
j = 1;

D1 = 0;
for p = pmin:pstep:pmax
    for k = 1:maxi
        for m = minm:maxm
            F = Fofi(maxi,k);
            CV = ConVoting(p,maxi,F);
            CVD = CappedVoteDeleg(p,maxi,F,m,cap);
            Delta = CV - CVD;
            if Delta > D1
                Result(j,1) = Delta;
                Result(j,2) = CV;
                Result(j,3) = CVD;
                Result(j,4) = k;
                Result(j,5) = maxi;
                Result(j,6) = m;
                Result(j,7) = p;
                j = j+1;
                D1 = Delta;
            end
            % for Mesh Plot
            %M(m,k) = m;
            %K(k,m) = k;
            %D(k,m) = Delta;
        end
    end
end

%M = M.';

T = table(Result(:,1),Result(:,2),Result(:,3),Result(:,4),Result(:,5),Result(:,6),Result(:,7));
T.Properties.VariableNames = {'Delta','CV','FVD','k','max i','m','p'};

clearvars -except T K M D


%% Plot
figure(1)
mesh(K,M,D)
xlabel('k')
ylabel('m')
zlh = zlabel('P(p) - P(p,m)');
%zlh.Position(1) = -16.6855;
%zlh.Position(2) = 109.8340;
%zlh.Position(3) = 0.0110;
xticks([0 10 20 30])
yticks([0 10 20 30 40 50])
zticks([0  0.025  0.050 ])
title('p = 0.6', 'Fontsize', 24)
zlim([0 0.05])
pbaspect([1 1 0.5])
view(-45,20)
set(gca,'fontsize', 18)


%%
function [F] = Fofi(maxi,k)
F = zeros(maxi,1);
F(k)=1;
end

%%
function [CV] = ConVoting(p,maxi,F)
CV = 0;
for i = 1:maxi
    if F(i) == 1
       value = 0;
       for k=0:i
           value = value + nchoosek(i,k) * p^k * (1-p)^(i-k) * gfunction(k,i-k);
       end
       CV = value;
    end
end
end

%%
function [FVD] = FreeVoteDeleg(p,maxi,F,m)
FVD = 0;
for i = 1:maxi
    if F(i) == 1
        value = 0;
        for k = 0:i
            value = value + nchoosek(i,k) * p^k * (1-p)^(i-k) * Gfunction(k,i-k,m);
        end
        FVD = value;
    end
end
end

%%
function [CVD] = CappedVoteDeleg(p,maxi,F,m,cap)
CVD = 0;
for i = 1:maxi
    if F(i) == 1
        value = 0;
        for k = 0:i
            value = value + nchoosek(i,k) * p^k * (1-p)^(i-k) * Gcapfunction(k,i-k,m,cap);
        end
        CVD = value;
    end
end

end

%%
function [G] = Gcapfunction(k,l,m,cap)
G = 0;
for h = 0:m
    G = G + nchoosek(m,h) * (k/(k+l))^h * (l/(k+l))^(m-h) * gfunction(min(k+h,cap*k),min(l+m-h,cap*l));
end
end

%%
function [G] = Gfunction(k,l,m)
G = 0;
for h = 0:m
    G = G + nchoosek(m,h) * (k/(k+l))^h * (l/(k+l))^(m-h) * gfunction(k+h,l+m-h);
end
end

%%
function [g] = gfunction(k,l)
if k>l
    g = 1;
elseif k==l
    g = 0.5;
else
    g = 0;
end
end

