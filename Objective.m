function z =Objective(x)
global ObjCost ObjSatisfaction Budget SatFactor CostFactor NEVAL NINALIDEVAL
NEVAL = NEVAL +1;
if  validateInput(x) == false
    z=[inf,inf]';
     NINALIDEVAL=NINALIDEVAL+1;
    return
end
% Objective 1 is Cost
z1 = sum(x.*ObjCost)/CostFactor;
z2 = (SatFactor-sum(x.*ObjSatisfaction))/SatFactor;
if z1*CostFactor > Budget
    z1 = inf;
    z2 = inf;
    NINALIDEVAL=NINALIDEVAL+1;
end
z = [z1,z2]';
end