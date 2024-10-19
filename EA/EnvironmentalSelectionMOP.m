function [Pop,FrontNo,SpCrowdDis] = EnvironmentalSelectionMOP(Pop,N,Cons,label,VAR)
% D=2;
M=2;

[FrontNo,MaxFNo] = NDSort(Pop,M,N,Cons,label,VAR);
Next = FrontNo < MaxFNo;
PopDec=unique(Pop,'rows','stable');
PopObj=Calobj(PopDec,label);
%% Calculate the special crowding distance in objective and decision space of each solution
SpCrowdDis_Obj = ModifiedCrowdingDistance(PopObj,FrontNo);
SpCrowdDis_Dec = ModifiedCrowdingDistance_dec(PopDec,FrontNo,PopObj);
avg_Obj = mean(SpCrowdDis_Obj);
avg_Dec = mean(SpCrowdDis_Dec);
can=avg_Dec/avg_Obj;
d=SpCrowdDis_Dec./SpCrowdDis_Obj;
index=find(d<=can);
a=length(index)/size(SpCrowdDis_Obj,2);
b=(size(SpCrowdDis_Obj,2)-length(index))/size(SpCrowdDis_Obj,2);
SpCrowdDis_Obj=mapminmax(SpCrowdDis_Obj);
SpCrowdDis_Dec=mapminmax(SpCrowdDis_Dec);
CrowdDis=a*SpCrowdDis_Dec+b*SpCrowdDis_Obj;
SpCrowdDis=mapminmax(CrowdDis);
% SpCrowdDis=max(SpCrowdDis_Obj,SpCrowdDis_Dec);
for i = 1 : MaxFNo
    Front   = find(FrontNo==i);
    Avg_Obj = mean(SpCrowdDis_Obj(Front));
    Avg_Dec = mean(SpCrowdDis_Dec(Front));
    replace = SpCrowdDis_Obj(Front)<=Avg_Obj & SpCrowdDis_Dec(Front)<=Avg_Dec;
    SpCrowdDis(Front(replace)) = min(SpCrowdDis_Obj(Front(replace)),SpCrowdDis_Dec(Front(replace)));
end
%% Select the solutions in the last front based on their crowding distances
Last     = find(FrontNo==MaxFNo);
[~,Rank] = sort(SpCrowdDis(Last),'descend');
Next(Last(Rank(1:N-sum(Next)))) = true;
%% Population for next generation
Pop = Pop(Next,:);
FrontNo    = FrontNo(Next)';
SpCrowdDis = SpCrowdDis(Next);
end
function CrowdDis = ModifiedCrowdingDistance(PopObj,FrontNo)
[N,M]    = size(PopObj);
CrowdDis = zeros(1,N);
Fronts   = setdiff(unique(FrontNo),inf);
for f = 1 : length(Fronts)
    Sol_in_Front = find(FrontNo==Fronts(f));
    Fmax  = max(PopObj(Sol_in_Front,:),[],1);
    Fmin  = min(PopObj(Sol_in_Front,:),[],1);
    for i = 1 : M    % The crowding distance is the cumulative  sum of the crowding distance of the two objective functions
        [~,Rank] = sortrows(PopObj(Sol_in_Front,i));
        CrowdDis(Sol_in_Front(Rank(1))) = CrowdDis(Sol_in_Front(Rank(1))) + 1;
        for j = 2 : length(Sol_in_Front)-1
            if Fmax(i) == Fmin(i)
                CrowdDis(Sol_in_Front(Rank(j))) = CrowdDis(Sol_in_Front(Rank(j)))+1;
            else
                CrowdDis(Sol_in_Front(Rank(j))) = CrowdDis(Sol_in_Front(Rank(j)))+(PopObj(Sol_in_Front(Rank(j+1)),i)-PopObj(Sol_in_Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
end
end
function CrowdDis = ModifiedCrowdingDistance_dec(Pop_dec,FrontNo,functionvalue)
N    = size(Pop_dec,1);
CrowdDis = zeros(1,N);
Fronts   = setdiff(unique(FrontNo),inf);
for f = 1 : length(Fronts)
    Sol_in_Front = find(FrontNo==Fronts(f));
    [~,Rank] = sortrows(functionvalue(Sol_in_Front,1));
    if length(Sol_in_Front)-1==0
        CrowdDis(Sol_in_Front(Rank(1))) = CrowdDis(Sol_in_Front(Rank(1)))+1;
    else
        CrowdDis(Sol_in_Front(Rank(1))) = CrowdDis(Sol_in_Front(Rank(1))) + 1;
        CrowdDis(Sol_in_Front(Rank(end))) = CrowdDis(Sol_in_Front(Rank(end))) + 1;
        for j = 2 : length(Sol_in_Front)-1
            CrowdDis(Sol_in_Front(Rank(j))) = (pdist2(Pop_dec(Sol_in_Front(Rank(j)),:),Pop_dec(Sol_in_Front(Rank(j-1)),:),'hamming')+ ...
                pdist2(Pop_dec(Sol_in_Front(Rank(j)),:),Pop_dec(Sol_in_Front(Rank(j+1)),:),'hamming'))/2;
        end
    end
end
end
