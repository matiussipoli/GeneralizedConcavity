clear;
        
%%%%%%%%
%Simulations Control
%%%%%%%%

%Simulation of the threshold (0 - Simulates, 1 - Takes pre defined for p=
%100, 2 - takes the narrowest case (r to the limit) (most low probability)
simTh = 1;

%value (for good starting point) (0 - Takes vcho, 1 - Takes the highest value by dichotomy)
validSol = 0; 
vcho = 2; %Value to take the appropriate treshold (only for Takes vcho)


%%%%%%%%
%Data Base
%%%%%%%%

%Specifying the problem variables
yhours = 8760;
tm = 10 ; % 10 time steps
vp_o = [100;99;95;90;86;94;92;80;86;78];
vp = vp_o(1:tm);

%ct = round(rand(1,tm)*20,2)';
ct_o = [13.11;0.710;16.98;18.68;13.57;15.15;14.86;7.84;13.11;3.42];
ct = ct_o(1:tm);

%Products

prn = 5;

pr1_o = ones(10,1)*2;
pr1 = pr1_o(1:tm);

pr2_o = [0;0;2.38;1.93;1.14;2.43;1.60;1.05;2.82;2.63];
pr2 = pr2_o(1:tm);

pr3_o = [0;0;0;0;0;2.15;0.920;4.52;4.90;2.19];
pr3 = pr3_o(1:tm);

pr4_o = [0;2.29;0;4.82;2.73;0;0;2.44;3.12;3.40];
pr4 = pr4_o(1:tm);

pr5_o = [3.4;0;0;2;0;0;0;3.57;4.51;4.45];
pr5 = pr5_o(1:tm);

%Consumption

dmMf = [100.25 104.93 106.23 107.69 108.16 110.60 115.17 118.70 119.02 122.62;%1
        99.220 103.02 105.45 107.34 109.05 108.82 106.36 101.14  97.43 102.84;
        100.29 100.36 102.74 106.18 107.61 109.09 112.96 115.92 117.04 120.21;%3
        100.36 98.860 100.94 102.51 103.93 105.32 109.73 110.98 112.09 116.70;
        101.35 104.75 105.89 109.64 113.45 119.13 121.60 124.29 128.39 131.00;%5
        101.75 100.98 101.15 102.65 104.24 107.81 112.79 112.28 114.71 116.93;
        101.67 103.81 106.38 110.12 114.21 115.76 118.34 119.36 122.46 128.03;%7
        102.81 104.92 104.91 108.00 106.30 112.54 112.71 113.58 117.35 118.80;
        102.58 108.06 108.55 110.37 114.19 115.86 114.36 117.54 120.21 122.61;%9
        103.13 105.08 101.28 110.05 115.39 107.79 116.52 113.82 107.18 96.31];

conCen = 3;
dm = dmMf(conCen,1:tm)';

%Price Dist

%Original data from pSeries

cd 'C:\Users\Guilherme\Google Drive\Doutorado\Estudos\20191001_Generalized Conc on Bilinear Strc\Matlab'

load('pSerie.mat');

qqplot(mean(pSerie'))

pMix = 120*ones(1,tm);

pSerieSize = size(pSerie,1);

pMixM = repmat(pMix,pSerieSize,1);

pMix_Se =  pMixM - pSerie(:,1:tm);

pmu = mean(pMix_Se);

pSig = cov(pMix_Se);

%estimating a and b

exp10 = ones(tm,1)*10;
afc = 2;
ap = exp10'*pmu'*afc;
bp = -exp10'*pmu'*afc;


tic;

pro = 0.85;
pr = 1;

%Coefs
muCoef = pmu';
SigCoef = pSig;
a = ap;
b = bp;


%Pre solution
alpLub = 1;
epsLub = (1-pro)/alpLub;

muLub = 0;
sigmaLub = 1;
pdLub = makedist('Normal','mu',muLub,'sigma',sigmaLub);

invEpsLub_a = icdf(pdLub,epsLub);
invEpsLub_b = icdf(pdLub,pro);
invEpsLub_haf = icdf(pdLub,epsLub/2);


SigCoefCho = chol(SigCoef);

cvx_begin quiet
    variables pn(tm) xp(tm) pr(prn) tLub;
    minimize( ct'*pn );
    subject to
        for t = 1 : tm
            vp(t) + pr(1)*pr1(t) + pr(2)*pr2(t) + pr(3)*pr3(t) + pr(4)*pr4(t) + pr(5)*pr5(t) + xp(t) == dm(t); 
            pn(t) - xp(t) >= 0;
            pn(t) >= 0;
        end

        for p = 1 : prn
            pr(p)  >= 0; 
        end

        %Restricoes CC - Lubin Two-sided linear chance constraints and
        %extensions 2016
        
        tLub >= norm(SigCoefCho*xp)
        ap - muCoef'*xp <= invEpsLub_a*tLub
        bp - muCoef'*xp >= invEpsLub_b*tLub
        ap - bp <= 2*invEpsLub_haf*tLub      
cvx_end

%evaluation functions
vx = sym('vx',[1 tm]);
vxt = vx';

%Originals
hat = (muCoef'*vxt - a)/sqrt(vxt'*SigCoef*vxt);
hbt = (b - muCoef'*vxt)/sqrt(vxt'*SigCoef*vxt);

gx = normcdf(hbt) + normcdf(hat) - 1;

gstart = eval(subs(gx,vxt',xp'));

fprintf('Time steps %d, Consumption scenario %d, Initial Probability %d, Solution cost %d, probability %d, time %d \n', tm, conCen, pro, ct'*(pn), gstart, toc);