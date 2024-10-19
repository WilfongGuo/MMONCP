function [PF,PS,Non_dominated_sol]=CMMONCT_EA(adjacent,N,MaxFE,experiment_num,label)
%%     input:
%                   adjacent            :    patient sample's PGIN
%                   N          :    population size
%                   MaxFE      :    the maximum number of function evaluation
%                   experiment_num  :    the number of algorithm rus
%       output:
%                   PF   :    the value of objective function of DNB obtained by MMPDNB
%                   PS   :    non dominated DNB obtained by MMPDNB


%% MDS
net=adjacent.subnetwork_adjacency;
Net=net+eye(size(net,2),size(net,2));  %% 网络&&约束
[z1,z2]=find(Net~=0);
z=[z2 z1];
N_1=length(Net);
A_adjacent=zeros(N_1);
for hz=1:size(z,1)

    A_adjacent(z(hz,2),z(hz,1))=1;
    A_adjacent(z(hz,1),z(hz,2))=1;

end
Cons=A_adjacent;
Cons(sum(Cons,2)==1,:)=[];
Cons=unique(Cons,'rows','stable');

%% DFVS
% A=adjacent.subnetwork_adjacency;
% [z2,z1]=find(A~=0);
% z=[z1 z2];
% N1=length(A);
% [N2,~]=size(z);
% %calculate the adjacency matrix of bipartite graph
% A_adjacent=zeros(N2,2*N1);
% for i=1:N2
%
%     A_adjacent(i,z(i,1))=1;
%     A_adjacent(i,z(i,2))=-1;
%     A_adjacent(i,N1+z(i,1))=N1;
%
% end
% Cons=A_adjacent;
% Cons(all(Cons==0,2),:)=[];%删除全零行
% Cons=unique(Cons,'rows');

%% NCUA
% Network=adjacent.subnetwork_adjacency;
% [z1,z2]=find(Network~=0);
% z=[z1 z2];
% N1=length(Network);
% [N2,~]=size(z);%calculate the adjacency matrix of bipartite graph
% A_adjacent=zeros(N2,N1);
% for i=1:N2
%
%     A_adjacent(i,z(i,1))=1;
%     A_adjacent(i,z(i,2))=1;
% end
% Cons=A_adjacent;
% Cons(all(Cons==0,2),:)=[];%删除全零行
% Cons=unique(Cons,'rows');



%% Set experiment parameters
Non_dominated_sol=cell(experiment_num,1);
N1=N*0.3;    % the population size of subpopulations
maxGen=ceil(MaxFE/N);
D=length(label); %维度
%%   EA

for EXP_NUM=1:experiment_num

    X=0;
    gen=1;
    %% Generate intial population.
    FE=0;   % the number of function evaluation
    Population1= randi([0,1],N,D);  %% 初始化
    Population2 = randi([0,1],N1,D);  %% 初始化
    Population3 = randi([0,1],N1,D);  %% 初始化

    f1=Calfunctionvalue(Population1,label);% Calculate the objective function of
    Zmin1       = min(f1,[],1);

    cons=Calcons(Population1,Cons);% Calculate the Constraint violation degree
    cons1=Calcons(Population2,Cons);
    cons2=Calcons(Population3,Cons);
    cons = [cons;cons1;cons2];cons(cons<0) = 0;VAR0 = max(sum(cons,2));
    if VAR0 == 0
        VAR0 = 1;
    end
    %%  评价
    [Population1, FrontNo, CrowdDis] = EnvironmentalSelectionMOP(Population1,N,Cons,label,VAR0);
    [~,D_Dec2]=CalFitness(Population2,label,Cons,VAR0,1);
    [~,D_Dec3]=CalFitness(Population3,label,Cons,VAR0,2);
    FE=FE+N+2*N1;

    %% Evolution
    while(FE<MaxFE)

        cp=(-log(VAR0)-6)/log(1-0.5);
        %             adjust the threshold
        if X < 0.5
            VAR=VAR0*(1-X)^cp;
        else
            VAR=0;
        end
        %%  产生子代
        MatingPool = [Population1(randsample(N,N),:)];
        [~,Mate2,Mate3]  = Neighbor4_Pairing_Strategy(MatingPool,Zmin1,1,label);
        rand_perm = randperm(N);
        Off1 = OperatorGAhalf([Mate2(rand_perm(1:N/4),:);Mate3(rand_perm(1:N/4),:)]);
        MatingPool = TournamentSelection_hamming(N/2,N,Population1,FrontNo,CrowdDis);
        Off2 = OperatorGAhalf(Population1(MatingPool,:));
        Offspring1=[Off1;Off2];

        ind2=TournamentSelection(2,N1,D_Dec2);
        ind3=TournamentSelection(2,N1,D_Dec3);
        Offspring2  = OperatorGAhalf(Population2(ind2,:));
        Offspring3  = OperatorGAhalf(Population3(ind3,:));

        FE=FE+size(Offspring1,1)+size(Offspring2,1)+size(Offspring3,1);
        hhhhh=Calobj(Offspring1,label);Zmin1       = min([Zmin1;hhhhh],[],1);

        %%   环境选择
        [Population1, FrontNo, CrowdDis] = EnvironmentalSelectionMOP(unique([Population1;Offspring1;Offspring2;Offspring3],'rows','stable'), N,Cons,label,VAR);
        [Population2,~,D_Dec2 ]  = EnvironmentalSelection2(unique([Population2;Offspring1;Offspring2;Offspring3],'rows','stable'),N1,label,Cons,1,VAR);
        [Population3,~,D_Dec3 ]  = EnvironmentalSelection2(unique([Population3;Offspring1;Offspring2;Offspring3],'rows','stable'),N1,label,Cons,2,VAR);

        X=X+1/maxGen;
        gen=gen+1;
    end
    %% Record non-dominated soutions
    f1=Calfunctionvalue(Population1,label);
    CV1=Calcons(Population1,Cons);
    ind=find(CV1==0);
    Population1=Population1(ind,:);
    objs=f1(ind,:);
    [FrontNo,~] = NDSort(objs,inf);
    outputpop=Population1(FrontNo==1,:,:);
    Non_dominated_sol{EXP_NUM}=outputpop;
end

Non_dominated_sol=cell2mat(Non_dominated_sol);

[pop,~,~]=unique(Non_dominated_sol,'rows');% De-duplication

functionvalue = Calfunctionvalue(pop,label);

[FrontNo,~] = M_non_domination_scd_sort(pop,functionvalue);

POP=pop(FrontNo==1,:);

FV=functionvalue(FrontNo==1,:);

[PF,mod_position]=sortrows(FV);

PS=POP(mod_position,:);
end

function CV=Calcons(pop,Cons)
ind=pop';
cv=Cons*ind;
CV=sum(cv==0,1)';
end

function f=Calfunctionvalue(pop,label)
%% 函数值
f=zeros(size(pop,1),2);


f(:,1)=sum(pop,2)/length(label);
f(:,2)=1-pop*label/length(label);
end

function index = TournamentSelection(K,N,varargin)
varargin    = cellfun(@(S)reshape(S,[],1),varargin,'UniformOutput',false);
[Fit,~,Loc] = unique([varargin{:}],'rows');
[~,rank]    = sortrows(Fit);
[~,rank]    = sort(rank);
Parents     = randi(length(varargin{1}),K,N);
[~,best]    = min(rank(Loc(Parents)),[],1);
index       = Parents(best+(0:N-1)*K);
end

function Offspring = OperatorGAhalf(Parent)
proM=1;
Parent1 = Parent(1:floor(end/2),:);
Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
[N,D]   = size(Parent1);
%% Genetic operators for binary encoding
% Uniform crossover
k = rand(N,D) < 0.5;

Offspring1    = Parent1;
Offspring2    = Parent2;
Offspring1(k) = Parent2(k);
Offspring2(k) = Parent1(k);
Offspring     = [Offspring1;Offspring2];
if size(Offspring,1)<size(Parent,1)
    ofd=Parent1(3,:);
    ofd(1,k(3,:))=Parent2(3,k(3,:));
    Offspring= [Offspring;ofd];
end
% Bit-flip mutation
Site = rand(size(Offspring,1),D) < proM/D;
Offspring(Site) = ~Offspring(Site);
end

function [Population,Fitness,D_Dec] = EnvironmentalSelection2(Population,N,label,Cons,HZ,VAR0)

%% Calculate the fitness of each solution
[Fitness,D_Dec ]=CalFitness(Population,label,Cons,VAR0,HZ);
%% Environmental selection
Next = zeros(size(Population,1),1).';
%     if sum(Next) < N
[~,Rank] = sort(Fitness);
Next(Rank(1:N)) = true;
% Population for next generation
Population = Population(Next==1,:);
Fitness=Fitness(Next==1);
D_Dec=D_Dec(Next==1);
[Fitness,rank] = sort(Fitness);
Population = Population(rank,:);
D_Dec=D_Dec(rank);
end

function [Fitness,D_Dec ]=CalFitness(Pop,label,Cons,VAR,hz)
objs=Calobj(Pop,label);
CV=Calcv(Pop,Cons);
PopCon   = sum(max(0,CV),2);
Feasible = PopCon <= VAR;
R = Feasible.*objs(:,hz) + ~Feasible.*(PopCon+1e10);
N = size(objs,1);
Distance = pdist2(Pop,Pop,"hamming");
Distance(logical(eye(length(Distance)))) = inf;
Distance = sort(Distance,2);
D_Dec = 1./(Distance(:,floor(sqrt(N)))+2);
Fitness = R;
end



function [FrontNo,MaxFNo] = NDSort(Pop,M,N,Cons,label,VAR)
nsort  = N;
PopObj=Calobj(Pop,label);
PopCon=Calcv(Pop,Cons);
Infeasible = PopCon>VAR;
PopObj(Infeasible,:) =repmat(max(PopObj,[],1),sum(Infeasible),1) +...
    repmat(sum(max(0,PopCon(Infeasible,:)),2),1,M);
[FrontNo,MaxFNo] = ENS_SS(PopObj,nsort);
end
function [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort)
[PopObj,~,Loc] = unique(PopObj,'rows');
Table   = hist(Loc,1:max(Loc));
[N,M]   = size(PopObj);
FrontNo = inf(1,N);
MaxFNo  = 0;
while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
    MaxFNo = MaxFNo + 1;
    for i = 1 : N
        if FrontNo(i) == inf
            Dominated = false;   %Dominated=0
            for j = i-1 : -1 : 1
                if FrontNo(j) == MaxFNo
                    m = 2;
                    while m <= M && PopObj(i,m) >= PopObj(j,m)
                        m = m + 1;
                    end
                    Dominated = m > M;
                    if Dominated || M == 2
                        break;
                    end
                end
            end
            if ~Dominated
                FrontNo(i) = MaxFNo;
            end
        end
    end
end
FrontNo = FrontNo(:,Loc);
end
