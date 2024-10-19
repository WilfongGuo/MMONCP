function [Mate1,Mate2,Mate3] = Neighbor4_Pairing_Strategy(MatingPop,Zmin,flag,label)
    objs=Calobj(MatingPop,label);
    Objs = objs;
    [Num,M] = size(Objs);
    Objs = (Objs - repmat(Zmin,Num,1));
    Objs = Objs./repmat(sqrt(sum(Objs.^2,2)),1,M);

    CosV = Objs * Objs';
    CosV = CosV - 3*eye(Num,Num);

    [~,SInd] = sort(-CosV,2);

    Nr = 10;
    Neighbor = SInd(:,1:Nr);

    Mate1 = MatingPop;

    P = ones(Num,2);
    for i = 1:Num
         P(i,1:2) = Neighbor(i,randsample(Nr,2));
    end

    Mate2 = MatingPop(P(:,1),:);
    Mate3 = MatingPop(P(:,2),:);
end

