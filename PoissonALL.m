%% Conventional Voting, Free Delegation, Capped Delegation with Poisson Distribution of n
clear
cap = 2;
p = 0.52;
MinM = 500; %>0
MaxM = 700;
StepM = 25;
MinN = 600;
MaxN = 800;
StepN = 25;

disp('start')
%{
    i=1;
for p=0.95:-0.05:0.55
    if p>=0.75
        MinM = 1;
        MaxM = 15;
        MinN = 1;
        MaxN = 15;
    elseif p>0.59
        MinM = 1;
        MaxM = 40;
        MinN = 1;
        MaxN = 40;
    else
        MinM = 90;
        MaxM = 100;
        MinN = 90;
        MaxN = 100;
    end
%}
%[FVD, CV, N, M, D, MaxDiff] = CalcALL(p, MinM, MaxM, StepM, MinN, MaxN, StepN);
[CVD, CV, N, M, D, MaxDiff] = CalcALLCap(p, MinM, MaxM, StepM, MinN, MaxN, StepN,cap);
%clearvars -except FVD CV N M D MaxDiff p CVD i
%Res(i,:)=MaxDiff;
%i=i+1;
%end

%% Search Max Diff with the Golden Section Algorithm
clear

p=0.55;
MinN = 100; %>2
MaxN = 180;
[M] = SearchMaxDelta(p,MinN,MaxN);
clearvars -except M

%% find max
for i = 1:Row
    mPos = find(D(i,:) == max(D(i,:)));
    mPos = mPos(1);
    mMaxDiff(i,1) = M(i,mPos);
    MaxDiff(i,1) = D(i,mPos);
end

Ex(:,1) = y';
Ex(:,2) = mMaxDiff;
Ex(:,3) = MaxDiff;

%% 3D Plot
figure(1)
mesh(N,M,D)
xlabel('n')
ylabel('m')
zlh = zlabel('P(n,p) - Pc(n,p,m)');
zlh.Position(1) = -16.2490;
zlh.Position(2) = 104.7951;
zlh.Position(3) = 0.0150;
xticks([0  20  40  60  80  100])
yticks([0  20  40  60  80  100])
zticks([0  0.01  0.02 0.03 ])
title(['p = ',num2str(p)], 'Fontsize', 24)
zlim([0 0.03])
pbaspect([1 1 0.7])
view(-45,20)
set(gca,'fontsize', 18)

%% Calculate CV, FVD and Diff
function[FVD, CV, N, M, D, MaxDiff] = CalcALL(p, MinM, MaxM, StepM, MinN, MaxN, StepN)
% This Function calculates all values for Conventional Voting, Free Vote
% Delegation and their Difference for m={MinM,...,MaxM} and
% n={MinN,..,MaxN} and p

%1. Calculate Boundaries of the Poisson Distribution for which all values will be calculated (max. k, l and min. k,l)
if p > 0.5
    klMax = poissinv(0.999999, MaxN*p);
    klMin = poissinv(0.000001, MinN*(1-p));
else
    klMax = poissinv(0.999999, MaxN*(1-p));
    klMin = poissinv(0.000001,MinN*p);
end

%2.Calculate G-Function
GF = PreCalc_G_FuncFree(klMin,klMax,MinM,MaxM);
disp('done with GF')

%3. Pre Calculate Poisson Values
% For p: Poissp(n,k)=poisspdf(n,k); Note that if we start with n or k >1
% all values will be zero in the matrix until n or k is reached
% As n and k can be 0 we have to start with Poissp(n+1 and k+1)
% Output: Poissp(n+1,k+1)=poisspdf(k,n*p)
for n=MinN:StepN:MaxN
    for k=klMin:klMax
        Poissp(n+1,k+1)=poisspdf(k,n*p);
    end
end
clear k
clear n

% For (1-p)
% For (1-p): Poiss1p(n,k)=poisspdf(n,k); Note that if we start with n or k >1
% all values will be zero in the matrix until n or k is reached
% As n and k can be 0 we have to start with Poissp(n+1 and k+1)
% Output: Poiss1p(n+1,k+1)=poisspdf(k,n*(1-p))
for n=MinN:StepN:MaxN
    for k=klMin:klMax
        Poiss1p(n+1,k+1)=poisspdf(k,n*(1-p)); 
    end
end
clear k
clear n

% Calculate all Values
b = 1;
for n = MinN:StepN:MaxN
    CVvalue = ConventionalVoting(n,klMin,klMax,Poissp,Poiss1p);
    t = 1;
    for m=MinM:StepM:MaxM
        FVDvalue = FreeVoteDelegation(n,klMin,klMax,Poissp,Poiss1p,GF,m);
        CV(b,t) = CVvalue;
        FVD(b,t) = FVDvalue;
        D(b,t) = CVvalue-FVDvalue;
        t = t+1;
    end
    b = b+1;
end


% Create Output Matrices
x = [MinM:StepM:MaxM];
y = [MinN:StepN:MaxN];
y = y';
%Number of rows
Row = floor((MaxN-MinN)/StepN)+1;
%Number of columns
Col = floor((MaxM-MinM)/StepM)+1;
% Create M Matrix - # Delgators
M = zeros(Row, Col);
for r = 1:Row
    M(r, :) = x;
end
% Create N Matrix - # Voters (expected)
N = zeros(Row, Col);
for col = 1:Col
    N(:,col) = y;
end


%find max Difference
maxPos = find(D == max(D(:)));
maxPos = maxPos(1);
MaxDiff(1,1) = D(maxPos);
MaxDiff(1,2) = CV(maxPos);
MaxDiff(1,3) = FVD(maxPos);
MaxDiff(1,4) = N(maxPos);
MaxDiff(1,5) = M(maxPos);
MaxDiff(1,6) = p;
end

%% Calculate CV, CVD and Diff
function[CVD, CV, N, M, D, MaxDiff] = CalcALLCap(p, MinM, MaxM, StepM, MinN, MaxN, StepN,cap)
% This Function calculates all values for Conventional Voting, Free Vote
% Delegation and their Difference for m={MinM,...,MaxM} and
% n={MinN,..,MaxN} and p

%1. Calculate Boundaries of the Poisson Distribution for which all values will be calculated (max. k, l and min. k,l)
if p > 0.5
    klMax = poissinv(0.999999, MaxN*p);
    klMin = poissinv(0.000001, MinN*(1-p));
else
    klMax = poissinv(0.999999, MaxN*(1-p));
    klMin = poissinv(0.000001,MinN*p);
end

%2.Calculate G-Function
GFcap = PreCalc_G_FuncCap(klMin,klMax,MinM,MaxM,cap);
disp('done with GF')

%3. Pre Calculate Poisson Values
% For p: Poissp(n,k)=poisspdf(n,k); Note that if we start with n or k >1
% all values will be zero in the matrix until n or k is reached
% As n and k can be 0 we have to start with Poissp(n+1 and k+1)
% Output: Poissp(n+1,k+1)=poisspdf(k,n*p)
for n=MinN:StepN:MaxN
    for k=klMin:klMax
        Poissp(n+1,k+1)=poisspdf(k,n*p);
    end
end
clear k
clear n

% For (1-p)
% For (1-p): Poiss1p(n,k)=poisspdf(n,k); Note that if we start with n or k >1
% all values will be zero in the matrix until n or k is reached
% As n and k can be 0 we have to start with Poissp(n+1 and k+1)
% Output: Poiss1p(n+1,k+1)=poisspdf(k,n*(1-p))
for n=MinN:StepN:MaxN
    for k=klMin:klMax
        Poiss1p(n+1,k+1)=poisspdf(k,n*(1-p)); 
    end
end
clear k
clear n

% Calculate all Values
b = 1;
for n = MinN:StepN:MaxN
    CVvalue = ConventionalVoting(n,klMin,klMax,Poissp,Poiss1p);
    t = 1;
    for m=MinM:StepM:MaxM
        CVDvalue = CappedVoteDelegation(n,klMin,klMax,Poissp,Poiss1p,GFcap,m);
        CV(b,t) = CVvalue;
        CVD(b,t) = CVDvalue;
        D(b,t) = CVvalue-CVDvalue;
        t = t+1;
    end
    b = b+1;
end


% Create Output Matrices
x = [MinM:StepM:MaxM];
y = [MinN:StepN:MaxN];
y = y';
%Number of rows
Row = floor((MaxN-MinN)/StepN)+1;
%Number of columns
Col = floor((MaxM-MinM)/StepM)+1;
% Create M Matrix - # Delgators
M = zeros(Row, Col);
for r = 1:Row
    M(r, :) = x;
end
% Create N Matrix - # Voters (expected)
N = zeros(Row, Col);
for col = 1:Col
    N(:,col) = y;
end


%find max Difference
maxPos = find(D == max(D(:)));
maxPos = maxPos(1);
MaxDiff(1,1) = D(maxPos);
MaxDiff(1,2) = CV(maxPos);
MaxDiff(1,3) = CVD(maxPos);
MaxDiff(1,4) = N(maxPos);
MaxDiff(1,5) = M(maxPos);
MaxDiff(1,6) = p;


end

%% Conventional Voting
function [CV] = ConventionalVoting(n,klMin,klMax,Poissp,Poiss1p) 
CV = 0;
for k=klMin:klMax
    for l=klMin:klMax
        CV = CV + Poissp(n+1,k+1) * Poiss1p(n+1,l+1) * gFunction(k,l);
        % Note: Poissp(n+1,k+1)=poisspdf(k,n*p), Poiss1p(n+1,k+1)=poisspdf(k,n*(1-p))
    end
end
clear k
clear l
end

%% Free Vote Delegation
function [FVD] = FreeVoteDelegation(n,klMin,klMax,Poissp,Poiss1p,GF,m) 
FVD = 0;
for k=klMin:klMax
    for l=klMin:klMax
        FVD = FVD + Poissp(n+1,k+1) * Poiss1p(n+1,l+1) * GF{1,m}(k+1,l+1);
        % Note: G{1,m}(k+1,l+1) = G-Function(k,l,m)
        % Note: Poissp(n+1,k+1)=poisspdf(k,n*p), Poiss1p(n+1,k+1)=poisspdf(k,n*(1-p))
    end
end
clear k
clear l
end

%% Capped Vote Delegation
function [CVD] = CappedVoteDelegation(n,klMin,klMax,Poissp,Poiss1p,GFcap,m) 
CVD = 0;
for k=klMin:klMax
    for l=klMin:klMax
        CVD = CVD + Poissp(n+1,k+1) * Poiss1p(n+1,l+1) * GFcap{1,m}(k+1,l+1);
        % Note: G{1,m}(k+1,l+1) = G-Function(k,l,m)
        % Note: Poissp(n+1,k+1)=poisspdf(k,n*p), Poiss1p(n+1,k+1)=poisspdf(k,n*(1-p))
    end
end
clear k
clear l
end

%% Pre calculate the G Function
function[G] = PreCalc_G_FuncFree(klMin,klMax,MinDel,MaxDel)
%Output: G{1,m}(k+1,l+1) = G-Function(k,l,m)

%Load Pascal Matrix to look up the Binomial Coefficient in the G Function
P = load('Pascal.mat');
P = P.P;

for m = MinDel:MaxDel
    for k = klMin:klMax
        for l = klMin:klMax
            % As k and l can be 0, we start with A(k+1 and l+1)
            A(k+1,l+1) = GFunctionFree(k,l,m,P);
        end
    end
    G{m} = A;
    clear A
end
clear k
clear l
clear m
end

%% G Function
function[G] = GFunctionFree(k, l, m,P)
G = 0;
res1 = 0;
A = k/(k+l);
B = l/(k+l);

%Calculate the first h for which g>0 can appear
F = floor((l+m-k)/2);
% If F<0, start with 0
if F < 0
    F = 0;
end

if k == l
    G = 0.5;
elseif k == 0
    G = 0;
elseif l == 0
    G = 1;
else
    for h = F:m
        if h > 300
            %For large m the Binomial part of the G function is approximated by
            %a normal distribution
            % check if approximation is reasonable
            if 1/(m^0.5) * abs((l/k)^0.5-(k/l)^0.5) >= 1/3
                disp('normal approx not appropriate')
            end
            res1 = normpdf(h,m*A,(m*A*(1-A))^0.5) * gFunction(k+h,l+m-h);
        else
            %For small m the G function is used directly
            res1 = P(m-h+1, h+1) * (A^h) * (B^(m-h)) * gFunction(k+h,l+m-h);
            % 
        end
        G = G + res1;
    end
end
clear k
clear l
clear h
clear res1
end

%% CAPPED Pre calculate the G Function
function[Gcap] = PreCalc_G_FuncCap(klMin,klMax,MinDel,MaxDel,cap)
%Output: G{1,m}(k+1,l+1) = G-Function(k,l,m)

%Load Pascal Matrix to look up the Binomial Coefficient in the G Function
P = load('Pascal.mat');
P = P.P;

for m = MinDel:MaxDel
    for k = klMin:klMax
        for l = klMin:klMax
            % As k and l can be 0, we start with A(k+1 and l+1)
            A(k+1,l+1) = GFunctionCap(k,l,m,P,cap);
        end
    end
    Gcap{m} = A;
    clear A
end
clear k
clear l
clear m
end

%% CAPPED G Function
function[Gcap] = GFunctionCap(k, l, m,P,cap)
Gcap = 0;
res1 = 0;
A = k/(k+l);
B = l/(k+l);

%Calculate the first h for which g>0 can appear
F = floor((l+m-k)/2);
% If F<0, start with 0
if F < 0
    F = 0;
end

if k == l
    Gcap = 0.5;
elseif k == 0
    Gcap = 0;
elseif l == 0
    Gcap = 1;
else
    for h = 0:m
        if h > 300
            %For large m the Binomial part of the G function is approximated by
            %a normal distribution
            % check if approximation is reasonable
            if 1/(m^0.5) * abs((l/k)^0.5-(k/l)^0.5) >= 1/3
                disp('normal approx not appropriate')
            end
            res1 = normpdf(h,m*A,(m*A*(1-A))^0.5) * gFunction(min(k+h,cap*k),min(l+m-h,cap*l));
        else
            %For small m the G function is used directly
            res1 = P(m-h+1, h+1) * (A^h) * (B^(m-h)) * gFunction(min(k+h,cap*k),min(l+m-h,cap*l));
            % 
        end
        Gcap = Gcap + res1;
    end
end
clear k
clear l
clear h
clear res1
end

%% g Function
function[g] = gFunction(k,l)
if k > l
    g = 1;
elseif k == l
    g = 0.5;
else
    g = 0;
end
end


%% Golden Section Search Algorithm
function[Max] = SearchMaxDelta(p,MinN,MaxN)
% The function looks for the maximal Difference between CV and FVD for
% n={MiN,..,MaxN} and for m = n-1

golden = (sqrt(5)-1)/2;

a = MinN;
b = MaxN;

x1 = round(a+(1-golden)*(b-a));
x2 = round(a+golden*(b-a));
[Free1, Conv1, N1, M1, D1] = CalcALL(p,x1-1,x1-1,1,x1,x1,1);  
[Free2, Conv2, N2, M2, D2] = CalcALL(p,x2-1,x2-1,1,x2,x2,1);
 
i = 1;
loop = 1;
while loop == 1
    if D1 > D2
        b = x2;
        x2=x1;
        x1 = floor(a+(1-golden)*(b-a));
        Max(i,1) = D1;
        Max(i,2) = Conv1;
        Max(i,3) = Free1;
        Max(i,4) = N1;
        Max(i,5) = M1;
        Max(i,6) = p;
        Dfinal = D1;
        Nfinal = N1;
        Mfinal = M1;
        [Free1, Conv1, N1, M1, D1] = CalcALL(p,x1-1,x1-1,1,x1,x1,1);  
        [Free2, Conv2, N2, M2, D2] = CalcALL(p,x2-1,x2-1,1,x2,x2,1);
    else
        a = x1;
        x1 = x2;
        x2 = ceil(a+golden*(b-a));
        Max(i,1) = D2;
        Max(i,2) = Conv2;
        Max(i,3) = Free2;
        Max(i,4) = N2;
        Max(i,5) = M2;
        Max(i,6) = p;
        Dfinal = D2;
        Nfinal = N2;
        Mfinal = M2;
        [Free1, Conv1, N1, M1, D1] = CalcALL(p,x1-1,x1-1,1,x1,x1,1);  
        [Free2, Conv2, N2, M2, D2] = CalcALL(p,x2-1,x2-1,1,x2,x2,1);
    end
    i = i+1;
    if abs(x1-x2)<=1
        loop = 0;
    end
end

% Because of the rounding the result might be not exactly the maximum,
% therefore we check for n-2, n-1, n+1, n+2
for k = -2:1:2
    new = Nfinal+k;
   [Freenew, Convnew, Nnew, Mnew, Dnew] = CalcALL(p,new-1,new-1,1,new,new,1);
   if Dnew>Dfinal
        Max(i,1) = Dnew;
        Max(i,2) = Convnew;
        Max(i,3) = Freenew;
        Max(i,4) = Nnew;
        Max(i,5) = Mnew;
        Max(i,6) = p;
        Dfinal = Dnew;
        Nfinal = Nnew;
        Mfinal = Mnew;
        Dfinal = Dnew;
   end
    
end


%look if found Delta is really the correct one in n
for j =-1:2:1
    newN = Nfinal+j;
    [Freed, Convd, N3, M3, D3] = CalcALL(p,Mfinal,Mfinal,1,newN,newN,1);
    if D3 > Dfinal
        disp('larger Diff found for m not equal to n-1');
    end
end

%look if found Delta is really the correct one in m
for j =-1:2:1
    newM = Mfinal+j;
    [Freed, Convd, Nd, Md, D3] = CalcALL(p,newM,newM,1,Nfinal,Nfinal,1);
    if D3 > Dfinal
        disp('larger Diff found for m not equal to n-1');
	end
end

end

%% Create the Pascal Matrix P
function[P] = CalcPascal()
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



