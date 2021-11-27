% Project Title: A Hybrid Multi-Objective Optimization Algorithm for Software Requirement Problem
% Developer: wathiq dukhan
% Contact Info: dralwathiq2012@gmail.com
%
function [NDS,HV,DeltaSpread] = HGABC
clear;

close all;
showFigure = true;
%% Problem Definition
global ObjCost P W ObjSatisfaction Budget ImplicationInteraction ...
    CompinationInteraction ExclusionInteraction InteractionEnabled ...
    SatFactor CostFactor RUN DIST NEVAL NINALIDEVAL DSTLIM TERMLIM ...
    pCrossover pMutation mu budgetRatio MaxIt nPop
NEVAL=0;
NINALIDEVAL=0;

% ObjCost=[1,4,2,3,4,7,10,2,1,3,2,5,8,2,1,4,10,4,8,4];
% P=[4,2,1,2,5,5,2,4,4,4,2,3,4,2,4,4,4,1,3,2;
%    4,4,2,2,4,5,1,4,4,5,2,3,2,4,4,2,3,2,3,1;
%    5,3,3,3,4,5,2,4,4,4,2,4,1,5,4,1,2,3,3,2;
%    4,5,2,3,3,4,2,4,2,3,5,2,3,2,4,3,5,4,3,2;
%    5,4,2,4,5,4,2,4,5,2,4,5,3,4,4,1,1,2,4,1]; 
% W=[1 4 2 3 4];
% 
% InteractionEnabled = true;
% ImplicationInteraction = [4 , 8;
%                           4 ,17;
%                           8 ,17;
%                           9 , 3;
%                           9 , 6;
%                           9 ,12;
%                           9 ,19;
%                           11,19];
% CompinationInteraction = [3 ,12;
%                           11,13;];
% ExclusionInteraction  = [];
%----------------------------
%dataset2
ObjCost=[16,19,16,7,19,15,8,10,6,18,15,12,16,20,9,4,16,2,9,3,...
         2,10,4,2,7,15,8,20,9,11,5,1,17,6,2,16,8,12,18,5,...
	     6,14,15,20,14,9,16,6,6,6,6,2,17,8,1,3,14,16,18,7,...
	     10,7,16,19,17,15,11,8,20,1,5,8,3,15,4,20,10,20,3,20,...
		 10,16,19,3,12,16,15,1,6,7,15,18,4,7,2,7,8,7,7,3 ];
P=[1,2,1,1,2,3,3,1,1,3,1,1,3,2,3,2,2,3,1,3,2,1,1,1,3,...
   3,3,3,1,2,2,3,2,1,2,2,1,3,3,2,2,2,3,1,1,1,2,2,3,3,...
   3,3,1,3,2,1,3,1,3,1,2,2,3,3,1,3,1,3,2,3,1,3,2,3,1,...
   1,2,3,3,1,2,1,3,1,2,2,2,1,3,2,2,3,1,1,1,2,1,3,1,1;
   3,2,1,2,1,2,1,2,2,1,2,3,3,2,1,3,2,3,3,1,3,3,3,2,3,...
   1,2,2,3,3,1,3,2,2,1,2,3,2,3,3,3,3,1,1,3,2,2,2,1,3,...
   3,3,1,2,2,3,3,2,1,1,1,3,2,3,1,2,1,2,3,1,1,3,1,3,2,...
   1,3,3,1,2,1,2,1,2,2,1,3,2,2,2,3,2,2,3,2,2,1,3,1,1;
   1,1,1,2,1,1,1,3,2,2,3,3,3,1,3,1,2,2,3,3,2,1,2,3,2,...
   3,3,1,3,3,3,2,1,2,2,1,1,3,1,2,1,3,1,3,3,3,3,1,3,2,...
   3,1,2,3,2,3,2,1,2,3,1,1,2,3,3,1,3,3,3,1,3,1,3,1,1,...
   2,3,3,1,2,1,2,3,2,3,1,2,2,3,3,3,3,2,1,1,2,3,3,2,3;
   3,2,2,1,3,1,3,2,3,2,3,2,1,3,2,3,2,1,3,3,1,1,1,2,3,...
   3,2,1,1,1,1,2,2,2,3,2,2,3,1,1,3,1,1,3,1,2,1,1,3,2,...
   2,1,3,2,1,3,3,1,2,3,2,2,3,3,3,1,2,1,2,1,2,3,3,2,2,...
   2,1,3,3,1,3,1,2,2,2,1,1,1,3,1,1,3,3,1,2,1,2,3,1,3;
   1,2,3,1,3,1,2,3,1,1,2,2,3,1,2,1,1,1,1,3,1,1,3,3,3,...
   2,2,3,2,3,1,1,3,3,2,2,1,1,2,1,3,1,1,2,1,2,3,3,2,2,...
   1,3,3,2,3,1,2,1,3,2,2,2,1,2,1,3,2,1,2,1,2,2,3,2,1,...
   3,2,3,1,3,3,2,1,2,2,2,2,1,3,3,3,1,1,3,1,3,3,3,3,3];
W=[1 5 3 3 1];
InteractionEnabled = true;
ImplicationInteraction = [2 , 24; 3 ,26; 3 ,27; 3 ,28; 3 , 29; 4 ,5; 6 ,7; 7,30; 10,32; 10,33; 14,32;
                          14 , 34; 14 ,37; 14 ,38; 16 ,39; 16 , 40; 17 ,43; 29 ,49; 29,50; 29,51; 30,52;
                        30 , 53; 31 ,55; 32 ,56; 32 ,57; 33 , 58; 36 ,61; 39 ,63; 40,64; 43,65; 46,68;
                          47 , 70; 55 ,79; 56 ,80; 57 ,80; 62 , 83; 62 ,84; 64 ,87];
                           
 CompinationInteraction = [21 ,22;
                           32 ,33;
                           46 ,47;
                           65 ,66];
ExclusionInteraction  = [];
CompinationInteraction= sort(CompinationInteraction,'descend');
ExclusionInteraction = sort(ExclusionInteraction,'descend');
ImplicationInteraction = sort(ImplicationInteraction,'descend');
%-----------------------------------
nVar=numel(ObjCost);             % Number of Decision Variables
elements = nVar;
DIST = ones(1,elements);
for i = 2 : DSTLIM
    e = round(elements/i);
    xe = ones(1,e);
    st =  round((elements -  e)/2);
    en =  st+e-1;
    DIST(st:en) = DIST(st:en)+ xe; 
end
DIST = DIST/sum(DIST);
ObjSatisfaction = W*P/sum(W);
SatFactor = sum(ObjSatisfaction);
CostFactor = sum(ObjCost);

Budget = sum(ObjCost)*budgetRatio;

CostFunction=@(x) Objective(x);      % Cost Function

%options = optimoptions(@gamultiobj,'PopulationType','bitstring');
%gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub,options);
%[x,fval] = gamultiobj(CostFunction,nVar);
VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax= 1;          % Upper Bound of Variables

% Number of Objective Functions
nObj=numel(CostFunction(randsample(VarMin:VarMax,nVar,true)));
L = round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)

a = 1;                    % Acceleration Coefficient Upper Bound

resultFile=['result/results-ratio-', num2str(budgetRatio)];
%% HGABC Parameters
%pCrossover=0.8;                         % Crossover Percentage
nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

%pMutation=0.2;                          % Mutation Percentage
nMutation=round(pMutation*nPop);        % Number of Mutants

%mu=0.01;                    % Mutation Rate

%sigma=0.1*(VarMax-VarMin);  % Mutation Step Size


%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];

pop=repmat(empty_individual,nPop,1);
% Initialize Best Solution Ever Found
BestSol.Cost = inf;%Budget
for i=1:nPop
       
    pop(i).Position=sample(nVar);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    if pop(i).Cost <= BestSol.Cost
        BestSol = pop(i);
    end
end
C = zeros(nPop, 1);
% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);
 
BestCost = zeros(MaxIt, 1);
%% HGABC Main Loop
F1=[];
for it=1:MaxIt
    
    newbee=repmat(empty_individual,nPop,1);
  % Recruited Bees
    for i = 1:nPop
        
        % Choose k randomly, not equal to i
        K = [1:i-1 i+1:nPop];
        k = K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi = a*unifrnd(-1, +1, VarSize);
        
        % New Bee Position
        newbee(k).Position = round((pop(i).Position+phi.*(pop(i).Position-pop(k).Position))/2);
       
        newbee(k).Cost=CostFunction(newbee(k).Position);
        newbee(k).Cost=CostFunction(newbee(k).Position);
        
        % Comparision
        if newbee(k).Cost <= pop(i).Cost
            pop(i) = newbee(k);
        else
            C(i) = C(i)+1;
        end
        
    end

        
    %Crossover
     popc=repmat(empty_individual,nPop/2,2);
    for k=1:nCrossover %insted by npop
     %  if randi([0,1])<pCrossover 
        i1=randi([1 nPop]);
        p1=pop(i1);
     %  else 
        i2=randi([1 nPop]);
        p2=pop(i2);
     %  end
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
        
        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
        
    end
    popc=popc(:);
    
    %   %%%%%%new old  % Mutation
    popmn=repmat(empty_individual,nMutation,1);
    for k=1:nMutation
        
        i=randi([1 nPop]);
        p=pop(i);
        
        popmn(k).Position=Mutate(p.Position,mu);
        
        popmn(k).Cost=CostFunction(popmn(k).Position);
      
    end
    popmn=popmn(:);
    %%%%%%end new old
%  Merge
    pop=[pop
          popmn
        popc]; %#ok
      %end add new wathiq
     
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    pop=SortPopulation(pop);

    newbee=repmat(empty_individual,nPop,1);
 
    % Onlooker Beess
    for i = 1:nPop
        
        % Choose k randomly, not equal to i
        K = [1:i-1 i+1:nPop];
        k = K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi = a*unifrnd(-1, +1, VarSize);
        
        % New Bee Position
        newbee(k).Position = round((pop(i).Position+phi.*(pop(i).Position-pop(k).Position))/2);
       
        newbee(k).Cost=CostFunction(newbee(k).Position);
        newbee(k).Cost=CostFunction(newbee(k).Position);
        
        % Comparision
        if newbee(k).Cost <= pop(i).Cost
            pop(i) = newbee(k);
        else
            C(i) = C(i)+1;
        end
        
    end
    
     %Crossover
     popc=repmat(empty_individual,nPop/2,2);
    for k=1:nCrossover %nCrossover insted by npop
     %  if randi([0,1])<pCrossover 
        i1=randi([1 nPop]);
        p1=pop(i1);
     %  else 
        i2=randi([1 nPop]);
        p2=pop(i2);
     %  end
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
        
        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
        
    end
     popc=popc(:);
     %-----------------------------
    popmn=repmat(empty_individual,nMutation,1);
    for k=1:nMutation
        
        i=randi([1 nPop]);
        p=pop(i);
        
        popmn(k).Position=Mutate(p.Position,mu);
        
        popmn(k).Cost=CostFunction(popmn(k).Position);
        
    end
    popmn=popmn(:);
 % Merge
    pop=[pop
         popc
         popmn]; 
     
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    pop=SortPopulation(pop);

    %scout bee phase
    L = round(0.6*nVar*nPop);
    popsc=repmat(empty_individual,L,1);
    for i = 1:3*nPop %L
%         if C(i)<= L = round(0.6*nVar*nPop);
        popsc(i).Position=sample(nVar);
         popsc(i).Cost=CostFunction(popsc(i).Position);
             C(i) = 0;
%         end
    end
    

     % Update Best Solution Ever Found
    for i = 1:nPop
        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end
      %end add new wathiq
     
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    pop=SortPopulation(pop);
    % REMOVE DUPLICATES
   
%     pop = pop(uniqueIds,:);
    costs = [pop.Cost]';
    nel = size(costs,1);
    uniqueIds=[];
    for ci = 1 : nel
        current = costs(ci,:);
        uniq = true;
        for j = ci+1:nel
            
            jth = costs(j,:);          
           cst = abs(current-jth);
           
           if all(cst<=TERMLIM)
               uniq = false;
           end
        end
        if uniq==true
            uniqueIds = [uniqueIds,ci];       
        end        
    end
    pop = pop(uniqueIds,:);
    
    needed = max(0,nPop - numel(pop));
    extra = repmat(empty_individual,needed,1);
    for vari=1:needed    
        extra(vari).Position=sample(nVar);
        extra(vari).Cost=CostFunction(extra(vari).Position);    
    end
    pop = [pop
        extra];
    
    % Truncate
    pop=pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F]=SortPopulation(pop);
    
    % Store F1
    F1=pop(F{1});
    
   
    % Plot F1 Costs
    if showFigure == true
           % figure(1);
        frame = PlotCosts(pop,F1);

        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        outfile = 'evolution.gif';

        % On the first loop, create the file. In subsequent loops, append.
        if it==1
            imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
        else
            imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
        end
        pause(0.01);
    end
    
end

%% Results
%figure(1);
%PlotCosts(pop,F1);
c = [F1.Cost]';

costs= c(:,1)*CostFactor;
statisfaction = (SatFactor-c(:,2)*SatFactor) * sum(W);
nF1 = numel(F1);
mnmxCost = minmax(costs');
mnmxSatsf=minmax(statisfaction');

positions = zeros(numel(F1),nVar);
for i=1:nF1 
   positions(i,:) =  F1(i).Position;
end
allcolumns = [costs,statisfaction,positions];

allUnique = unique(allcolumns,'rows');
DeltaSpread = spread( unique(c,'rows'));
cunique = abs(unique(c,'rows'));
%if(cunique(1,1)==0.0)
% cunique = cunique(2:end,:);
%end
%cunique(:,2) = 10 ./ abs(cunique(:,2));
rHyper = max(cunique);
%cunique(:,2)= rHyper(2)-cunique(:,2);
%rHyper = max(cunique);
%Hypervolume indicator
HV = hypervolume(cunique,rHyper,100000);%100 or 10 defult 1000 ???? ??? ?????? ??? ???????? 
%HVR = HV* rHyper(1)*rHyper(2);
nuniq = numel(allUnique(:,1));
NDS = nuniq;
summary = ['min-cost,',num2str(mnmxCost(1)),',min-satisfactoin,',num2str(mnmxSatsf(1)), ...
           '\r\nmax-cost,',num2str(mnmxCost(2)),',max-satisfactoin,',num2str(mnmxSatsf(2)),...
		   ',HV,' ,num2str(HV),',Spread,' num2str(DeltaSpread),',NDS,' num2str(NDS)];
rs = sprintf('r%.0f,' , 1:nVar);
rs = rs(1:end-1);
data = [summary '\r\nCost,Satisfaction,' rs];
for i = 1: nuniq
    row  = sprintf('%.0f,' , allUnique(i,:));
    row = row(1:end-1);
    data = [data,num2str('\n'),row];
end
sprintf('%.0f,' , allUnique(i,:));
fileID = fopen([resultFile ,'-',  num2str(RUN) , '.csv'],'w');
fprintf(fileID,data);
fclose(fileID);
%changed by wathiq 2021
%csvwrite(resultFile,data)
% rsultsheet = 1;
% xlswrite(resultFile,allUnique,rsultsheet,'B5')
% xlswrite(resultFile,{'Budget',budgetRatio},rsultsheet,'A1')
% xlswrite(resultFile,{'min',mnmxCost(1),mnmxSatsf(1)},rsultsheet,'A2')
% xlswrite(resultFile,{'max',mnmxCost(2),mnmxSatsf(2)},rsultsheet,'A3')
% xlswrite(resultFile,{'Cost','Satisfaction'},rsultsheet,'B4')
% inputsheet = 2;

% NDS
% HV
% DeltaSpread
% 
% 
% NEVAL 
% NINALIDEVAL
end




