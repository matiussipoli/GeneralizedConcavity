clear;
    
%%%%%%%%
%Simulations Control
%%%%%%%%

nameSim='simTm10Con10p35v1';

%Simulation of the threshold (0 - Simulates, 1 - Takes pre defined for p=
%100 from file
simTh = 1;

%For the min poin : minPoin (0 - takes vmin, 1 - takes the narrowest case (r to the limit) (most low probability)
minPoint = 1;
%vmin = 200; %only in case of minPoin = 0

%For the slater point = value (for good starting point) (0 - Takes vcho, 1 - Takes the highest value by dichotomy->)
%0 - Takes vcho is used when the solution is ensured
%1 - Is used when one is looking for a slater point by a algorithm
validSol = 1; 
vcho = 2; %Value to take the appropriate treshold (only for 0 - Takes vcho)
%vcho = 100; %Value to take the appropriate treshold (only for 0 - Takes vcho)

% rcho -> define the r for the TR algorithm
rcho = -100;

%%%%%%%%
%Data Base
%%%%%%%%

tm = 10 ; % 10 time steps
conCen = 10; % demand pattern
pro = 0.35; % probability

%Specifying the problem variables
yhours = 8760;
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

%conCen = 10;
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

%pro = 0.85;

pr = 1;

%Coefs
muCoef = pmu';
SigCoef = pSig;
a = ap;
b = bp;

%%%%%%
%Functions Definitions
%%%%%%

%src = sqrt(muCoef'*inv(pSig)*muCoef); % right side of the initial
%inequality  <- It was with this sqrt before, I don't know why
src = muCoef'*inv(pSig)*muCoef; % right side of the initial inequality  


if simTh == 0

    % End of definition of functions
    syms  r t p;

    % Defining limit theta

    sRtpa = symfun((r-1)*exp(-1/2*t)/(sqrt(2*pi)*p),[t r p]);

    deltRtpa =   symfun(t - sqrt(t)*sRtpa(t,r,p) - 1,[t r p]);

    deltComp = symfun(t*(deltRtpa(t,r,p)-2)/(deltRtpa(t,r,p)+2),[t r p]);
    
    pb = 1; %considered probability for the dichotomy process
    rini = 1; % initial r
    stIni = 0.1; % step of the search
    sr = 200;

    thetVec = zeros(sr,1);
    rVec = zeros(sr,1);

    %Troca de sinal

    %Achar pontos extremos
    for cSt = 1:sr

        rbus = rini - cSt + 1;

        dPla = symfun(deltComp(t,rbus,pb),t);

        flag = 0;

        t1 = 0;
        t2 = 0;

        while flag == 0

            t1 = t2;
            t2 = t2 + stIni;
            if( dPla(t2) > src)
                flag =1;
            end

        end

        %Start dicotomy

        flag = 0;
        lim = 0.0001;

        while flag == 0

            tmid = (t1 + t2)/2;

            if(dPla(tmid) > src)
                t2 = tmid;
            else
                t1 = tmid;
            end

            if (abs(dPla(t2) - src)/src < lim)
                flag = 1;
            end
        end

        thetVec(cSt) = t2 ;
        rVec(cSt) =  rbus;

    end
    
    %tha = thetVec(vcho);
    %thb = thetVec(vcho);
    %rcho = rVec(vcho);
     
elseif simTh == 1
    
    load('thetVec100.mat');
    load('rVec100.mat');
    
    rVec = rVec100;
    thetVec = thetVec100;   
        
end


if minPoint == 0
    
    tam = size(thetVec);
    tha = thetVec(tam(1));
    thb = thetVec(tam(1));

else
    
    tha = src;
    thb = src;

end


pSig = cov(pMix_Se);


%%%%%%%%%%%
%Defining the feasible set 
%%%%%%%%%%%

tha = src;
thb = src;

pSiga = pSig*tha;
pSigb = pSig*thb;

pSigS = sqrtm(pSig);
pSigSa = sqrtm(pSiga);
pSigSb = sqrtm(pSigb);

cvx_begin quiet
    variables pn(tm) xp(tm) pr(prn);
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

        %Restricoes CC
        { pSigSa*xp, muCoef'*xp - ap } <In> lorentz(tm);
        { pSigSb*xp, bp - muCoef'*xp } <In> lorentz(tm);

cvx_end


vx = sym('vx',[1 tm]);
vxt = vx';

%Originals
hat = (muCoef'*vxt - a)/sqrt(vxt'*SigCoef*vxt);
hbt = (b - muCoef'*vxt)/sqrt(vxt'*SigCoef*vxt);
gx = normcdf(hbt) + normcdf(hat) - 1;

pmin = eval(subs(gx,vxt',xp'));
fprintf('Minimal Probability: %d \n',pmin);


%%%%%%%%%%%
%Defining the miniminal set by simulation 
%%%%%%%%%%%

tam = size(thetVec);
tha = thetVec(tam(1));
thb = thetVec(tam(1));

pSiga = pSig*tha;
pSigb = pSig*thb;

pSigS = sqrtm(pSig);
pSigSa = sqrtm(pSiga);
pSigSb = sqrtm(pSigb);

cvx_begin quiet
    variables pn(tm) xp(tm) pr(prn);
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

        %Restricoes CC
        { pSigSa*xp, muCoef'*xp - ap } <In> lorentz(tm);
        { pSigSb*xp, bp - muCoef'*xp } <In> lorentz(tm);

cvx_end


vx = sym('vx',[1 tm]);
vxt = vx';

%Originals
hat = (muCoef'*vxt - a)/sqrt(vxt'*SigCoef*vxt);
hbt = (b - muCoef'*vxt)/sqrt(vxt'*SigCoef*vxt);
gx = normcdf(hbt) + normcdf(hat) - 1;

pminSim = eval(subs(gx,vxt',xp'));
fprintf('Minimal Probability Simulation: %d \n',pminSim);


%%%%%%%%%%%
%Defining probability for r=0
%%%%%%%%%%%

r0pos = find(rVec==0);

tha = thetVec(r0pos);
thb = thetVec(r0pos);

pSig = cov(pMix_Se);

pSiga = pSig*tha;
pSigb = pSig*thb;

pSigS = sqrtm(pSig);
pSigSa = sqrtm(pSiga);
pSigSb = sqrtm(pSigb);

cvx_begin quiet
    variables pn(tm) xp(tm) pr(prn);
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

        %Restricoes CC
        { pSigSa*xp, muCoef'*xp - ap } <In> lorentz(tm);
        { pSigSb*xp, bp - muCoef'*xp } <In> lorentz(tm);

cvx_end


vx = sym('vx',[1 tm]);
vxt = vx';

%Originals
hat = (muCoef'*vxt - a)/sqrt(vxt'*SigCoef*vxt);
hbt = (b - muCoef'*vxt)/sqrt(vxt'*SigCoef*vxt);
gx = normcdf(hbt) + normcdf(hat) - 1;

plog = eval(subs(gx,vxt',xp'));
fprintf('Probability of plog: %d \n',plog);

%Defining Slater Solution

itMax = 0;
itDic = 0;
dicFase = 1;
it1 = 0;
it2 = 0;
dicEps = 0.00001;

while itMax == 0

% General variables

if validSol ==  0
    % Just one calculation for the defined vcho
    tha = thetVec(vcho);
    thb = thetVec(vcho);
    
    itMax = 1; 
else
    % Dicotomy procedure
    
    itDic = itDic + 1;
    
    if dicFase == 1
       
        tha = thetVec(itDic);
        thb = thetVec(itDic);

    else 
        itm = (it1 + it2)/2;
        tha = itm;
        thb = itm;
    end
        
end

    pSiga = pSig*tha;
    pSigb = pSig*thb;
    
    pSigS = sqrtm(pSig);
    pSigSa = sqrtm(pSiga);
    pSigSb = sqrtm(pSigb);

    cvx_begin quiet
        variables pn(tm) xp(tm) pr(prn);
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

            %Restricoes CC
            { pSigSa*xp, muCoef'*xp - ap } <In> lorentz(tm);
            { pSigSb*xp, bp - muCoef'*xp } <In> lorentz(tm);

    cvx_end
    
    if strcmp(cvx_status,'Infeasible') && dicFase == 1
       
        it1 = it2;
        it2 = tha;

    elseif ~strcmp(cvx_status,'Infeasible') && dicFase == 1;
       
        it1 = it2;
        it2 = tha;
        
        dicFase = 2;
        
        if itDic == 1
            itMax =1;
        end
    
    elseif strcmp(cvx_status,'Infeasible') && dicFase == 2;
        
        it1 = itm;
        
    elseif ~strcmp(cvx_status,'Infeasible') && dicFase == 2;
        
        it2 = itm;
        
        if abs(it2 - it1)/it2 < dicEps
       
            itMax =1;
       
        end
        
    end
    
    
end


% Starting the algorithm resolution - chosen method, barrier method with
% TR method.
tic;

%Defining the functions and Sets
vx = sym('vx',[1 tm]);
vxt = vx';

%Originals
hat = (muCoef'*vxt - a)/sqrt(vxt'*SigCoef*vxt);
hbt = (b - muCoef'*vxt)/sqrt(vxt'*SigCoef*vxt);
dhat = -1/sqrt(vxt'*SigCoef*vxt)*(-muCoef + ((muCoef'*vxt - a)*SigCoef*vxt)/(vxt'*SigCoef*vxt));
dhbt = -1/sqrt(vxt'*SigCoef*vxt)*(muCoef + ((b - muCoef'*vxt)*SigCoef*vxt)/(vxt'*SigCoef*vxt));
ddthat = (1/(vxt'*SigCoef*vxt))*(muCoef*muCoef' - ((muCoef'*vxt - a)/(vxt'*SigCoef*vxt)*(muCoef*vxt'*SigCoef + SigCoef*vxt*muCoef'))+((muCoef'*vxt - a)/(vxt'*SigCoef*vxt))^2*SigCoef*vxt*vxt'*SigCoef);  
ddthbt = (1/(vxt'*SigCoef*vxt))*(muCoef*muCoef' - ((b - muCoef'*vxt)/(vxt'*SigCoef*vxt)*(muCoef*vxt'*SigCoef + SigCoef*vxt*muCoef'))+((b - muCoef'*vxt)/(vxt'*SigCoef*vxt))^2*SigCoef*vxt*vxt'*SigCoef);  
d2that = (1/(vxt'*SigCoef*vxt))^(3/2)*(- SigCoef*vxt*muCoef' - muCoef*vxt'*SigCoef - (muCoef'*vxt - a)*SigCoef + 3*((muCoef'*vxt - a)/(vxt'*SigCoef*vxt))*SigCoef*vxt*vxt'*SigCoef);
d2thbt = (1/(vxt'*SigCoef*vxt))^(3/2)*(- SigCoef*vxt*muCoef' - muCoef*vxt'*SigCoef - (b - muCoef'*vxt)*SigCoef + 3*((b - muCoef'*vxt)/(vxt'*SigCoef*vxt))*SigCoef*vxt*vxt'*SigCoef);

%main functions

r = rcho;
%r = -98;
% I do not remember the right r that was fixed here. maybe -98 (miniminum
% that we calculate
gx = normcdf(hbt) + normcdf(hat) - 1;
if r ~= 0 
    fx = gx^r;
else
    fx = log(gx);
end

kx = fx - (pro)^r;
dtkx = r*(gx)^(r-1)*(exp(-1/2*hat^2)/sqrt(2*pi)*dhat + exp(-1/2*hbt^2)/sqrt(2*pi)*dhbt);
d2fax = r*(gx)^(r-1)*(exp(-1/2*hat^2)/sqrt(2*pi))*(((r-1)*exp(-1/2*hat^2)/(sqrt(2*pi)*gx))*dhat*dhat' -hat*dhat*dhat' + d2that + ((r-1)*exp(-1/2*hbt^2)/(sqrt(2*pi)*gx))*dhbt*dhat');
d2fbx = r*(gx)^(r-1)*(exp(-1/2*hbt^2)/sqrt(2*pi))*(((r-1)*exp(-1/2*hbt^2)/(sqrt(2*pi)*gx))*dhbt*dhbt' -hbt*dhbt*dhbt' + d2thbt + ((r-1)*exp(-1/2*hat^2)/(sqrt(2*pi)*gx))*dhat*dhbt');
d2kx = d2fax + d2fbx;

gstart = eval(subs(gx,vxt',xp'));
fprintf('Initial Solution Probability: %d \n',gstart);

if gstart < pro || strcmp(cvx_status,'Infeasible')
    fprintf('ALERT: Invalid initial point \n');
    return;
end

fprintf('Target probability: %d \n',pro);

%compositions
phi = -log(-kx);
dphi = 1/(-kx)*dtkx;
d2phi =  1/(kx^2)*dtkx*dtkx' + 1/(-kx)*d2kx;

vphi = eval(subs(phi,vxt',xp'));

%%%%%%%%%%
%Starting the barrier algorithm
%%%%%%%%%%

ft = 1; %first iteration
cvgence = zeros(100,9);
itB = 0;
itTR = 0;
itTot = 0;

%initial values for barrier
tpb = 1;
mu = 1.5;
epsb = 0.01;
mp = 1;

%TR algorithm
deltaIni = 5;
ro_0 = 0.5;
c0 = 1;
ver = 1;
epsn = 0.25;

%Pre verification
verProb = eval(subs(gx,vxt',xp'));
cvgence(1,7) = ct'*pn;
cvgence(1,8) = verProb;

%Starting the barrier algorithm

while mp/tpb > epsb

    itB = itB + 1;
    
    lambx = 1;
    
    ft = 1; %first iteration
    
    delta = deltaIni;
    
    ver = 1;
    
    %Starting the Newton's method
    
    while ver > epsn
    
        itTot = itTot + 1;
        itTR = itTR + 1;
    
        cvgence(itTot+1,1) = itTot;
        cvgence(itTot+1,2) = itTR;
        cvgence(itTot+1,3) = itB;
        
        fbal = zeros(tm,1);
        feq = zeros(tm,1);

        %Setting the fixed variables

        if ft == 1
            vxpb = xp;
            vpnb = pn;
            vprb = pr;
            ft = 0;
        end
            
        for t = 1 : tm
            fbal(t) = vxpb(t);
            feq(t) = - vxpb(t);      
        end


        %Fixing functions
        vphi = eval(subs(phi,vxt',vxpb'));
        vdphi = eval(subs(dphi,vxt',vxpb'));
        vd2phi = eval(subs(d2phi,vxt',vxpb'));

        cvx_begin quiet
            variables pnb(tm) xpb(tm) prb(prn);
            minimize (tpb*(ct'*(pnb)) + vphi + vdphi'*xpb + 1/2*xpb'*vd2phi*xpb);
            subject to
                for t = 1 : tm
                    vp(t) + prb(1)*pr1(t) + prb(2)*pr2(t) + prb(3)*pr3(t) + prb(4)*pr4(t) + prb(5)*pr5(t) + xpb(t) + fbal(t) == dm(t); 
                    pnb(t) - xpb(t) + feq(t) >= 0;
                    pnb(t) >= 0;
                end

                for p = 1 : prn
                    prb(p)  >= 0; 
                end

                norm(xpb,2) <= delta;
                
        cvx_end
       
        ared = eval(subs(phi,vxt',vxpb')) - eval(subs(phi,vxt',vxpb' + xpb'));
        pred = eval(subs(phi,vxt',vxpb')) - (vphi + vdphi'*xpb + 1/2*xpb'*vd2phi*xpb);
        
        rtest = ared/pred;

        %testing the probability
        gxtest = eval(subs(gx,vxt',vxpb' + xpb'));
        
        %if (rtest <= c0) || imag(rtest) ~= 0 || (gxtest <= pmin)  % if predidiction is too big, we must reduce the step
        if (rtest <= c0) || imag(rtest) ~= 0 % if predidiction is too big, we must reduce the step
            
            vpnb = vpnb;
            vxpb = vxpb;
            vprb = vprb;
            
            delta = delta*ro_0;
            
        else
            
            vpnb = pnb;
            vxpb = vxpb + xpb;
            vprb = prb;
            
        end 
        
        %testing the probability
        %gxfinal = eval(subs(gx,vxt',vxpb' + xpb'));
        gxfinal = eval(subs(gx,vxt',vxpb'));
        
        ver = norm(xpb,2);
        cvgence(itTot+1,4) = norm(xpb,2);
        cvgence(itTot+1,5) = tpb;
        cvgence(itTot+1,6) = cvx_optval;
        cvgence(itTot+1,7) = ct'*(pnb);
        cvgence(itTot+1,8) = gxfinal;
        cvgence(itTot+1,9) = toc;
       
        fprintf('Solution cost %d, probability %d, iterations %d, time %d \n', ct'*(pnb), gxfinal, itTot,toc);
     
    end % end of newthon's method
        
    %barrier method settings
    tpb = mu*tpb;
    
    xp = vxpb;
    pn = vpnb;
    pr = vprb;
    
    itTR = 0;

end % end of barrier method

fprintf('Solution cost %d, probability %d, iterations %d, time %d \n', ct'*(pnb), gxfinal, itTot,toc);

% printing solutions

inputsSim = {'tm', num2str(tm); 'Demand', num2str(conCen); 'Probability', num2str(pro)}; 
col_header = {'itTot', 'itTR', 'itBarrier', 'norm_xpb', 'tpb', 'cvx_optval', 'ctXpnb', 'cvx_optval', 'time_S'}; 
FileDestination = strcat(pwd,'\Saidas\',nameSim,'.xls');

delete(FileDestination)
xlswrite(FileDestination,{'Minimal Probability -Inf: ', num2str(pmin)},'Sheet1','L1'); 
xlswrite(FileDestination,{'Minimal Probability Sim: ', num2str(pminSim)},'Sheet1','L2');
xlswrite(FileDestination,{'Probability of plog: ', num2str(plog)},'Sheet1','L3');
xlswrite(FileDestination,{'Initial Solution Probability: ', num2str(gstart)},'Sheet1','L4');
xlswrite(FileDestination,{'Target probability: ', num2str(pro)},'Sheet1','L5');
xlswrite(FileDestination,inputsSim,'Sheet1','O1');
xlswrite(FileDestination,col_header,'Sheet1','A1'); 
xlswrite(FileDestination,cvgence,'Sheet1','A2');
