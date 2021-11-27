% Project Title: Project Title: A Hybrid Multi-Objective Optimization Algorithm for Software Requirement Problem
function main
global RUN DSTLIM TERMLIM pCrossover pMutation mu budgetRatio MaxIt nPop
    function dualfprintf(fid,message)
        fprintf(fid,[message '\n']);
        fprintf([message '\n']);
    end
rng(100)
runs =30;

budgetRatio = 1;%1; %0.3%0.5%0.7;
MaxIt = 30;     % Maximum Number of Iterations
nPop  = 100;      % Population Size
TERMLIM = -1; % if you want to filter the duplicate make this 0
DSTLIM = 1; % 1 if you want to have more solutios in the middle of the graph, increase this value, say 4
pCrossover = 0.80 ;%
pMutation =0.20 ;%
mu =0.01; %

fileID = fopen('main_non_reg_out.txt','w');
%termlim = [0,0.0001,0.0005,0.001,0.005,0.01]
NDS = zeros(1,runs);
HV = zeros(1,runs);
DeltaSpread = zeros(1,runs);
for i = 1: runs
    RUN = i; % pCrossover : 0.8, pMutation : 0.5 , mu : 0.01
    dualfprintf(fileID,['RUN : ' num2str(i) ' , TERMLIM : ' num2str(TERMLIM) ' , DSTLIM : ' num2str(DSTLIM), ', pCrossover : ' num2str(pCrossover) ', pMutation : ' num2str(pMutation) ' , mu : ' num2str(mu)]);
    tic
    [NDS(i),HV(i),DeltaSpread(i, IGD)] = HGABC;
    sec = toc;
    dualfprintf(fileID,['IGD : ' num2str(IGD(i)) ', NDS : ' num2str(NDS(i)) ' , HV : ' num2str(HV(i)) ' , SPREAD : ' num2str(DeltaSpread(i)),' RT: ' num2str(sec)]);
    dualfprintf(fileID,'----------------------------------------------------------------');
    
end

dualfprintf(fileID,['NDS : ' num2str(NDS)])
NDSavg = mean(NDS);
NDSstd = std(NDS);
dualfprintf(fileID,['NDSavg : ' num2str(NDSavg)])
dualfprintf(fileID,['NDSstd : ' num2str(NDSstd)])
HVavg = mean(HV);
HVstd = std(HV);
dualfprintf(fileID,['HV : ' num2str(HV)])
dualfprintf(fileID,['HVavg : ' num2str(HVavg)])
dualfprintf(fileID,['HVstd : ' num2str(HVstd)])
DeltaSpreadavg = mean(DeltaSpread);
DeltaSpreadstd = std(DeltaSpread);
dualfprintf(fileID,['DeltaSpread : ' num2str(DeltaSpread)])
dualfprintf(fileID,['DeltaSpreadavg : ' num2str(DeltaSpreadavg)])
dualfprintf(fileID,['DeltaSpreadstd : ' num2str(DeltaSpreadstd)])
fclose(fileID);
end

