function index = TournamentSelection_hamming(K,N,Population,FrontNo,SpCrowdDis)
index=zeros(K,1);
hamming_dist=pdist2(Population,Population,"hamming");
[~,site2]=sort(hamming_dist,2);
K1=randperm(N,K);  % parent1 number
K2=site2(K1,2);    % parent2 number
for i=1:K
    if FrontNo(K1(i))<FrontNo(K2(i))
        index(i)=K1(i);
    elseif FrontNo(K1(i))==FrontNo(K2(i))
        if SpCrowdDis(K1(i))>=SpCrowdDis(K2(i))
            index(i)=K1(i);
        else
            index(i)=K2(i);
        end
    else
        index(i)=K2(i);
    end
end
end
