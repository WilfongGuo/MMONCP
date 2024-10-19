function [subnetwork_genes,subnetwork_adjacency] = construct_PGIN(i,Normal, Tumor,gene_list,Net,Ref)
     
    %for the i-th sample
    %construct the tumor SSN
    
    sample_tumor=Tumor(:,i);
    [R0,P]=SSN(sample_tumor,Ref);
    
    P(isnan(P))=0;
    P(abs(P)>=0.05)=0;
    P(P~=0)=1;

    %construct the normal SSN
    clear  sample_tumor 
    sample_normal=Normal(:,i);
    [R1,P1]=SSN(sample_normal,Ref);
    clear  sample_normal 
    
    P1(isnan(P1))=0;
    P1(abs(P1)>=0.05)=0;
    P1(P1~=0)=1;
    
    C=P-P1; 
    R=(log2(abs(R0./R1)));
    R(isnan(R))=0;
    
    %sample_network{i,1}=C;
    D=abs(C).*R;
   %*****************save network data**************************8

   CC=D.*Net;
   CC(isnan(CC))=0;CC(CC==inf)=0;

   k=sum(CC);%
   [~,a]=find(k~=0);
    
   subnetwork_nodes=a';
   subnetwork_genes=gene_list(subnetwork_nodes);

   scores=k(a)';

   subnetwork_adjacency0=CC(subnetwork_nodes,:);
   subnetwork_adjacency=subnetwork_adjacency0(:,subnetwork_nodes);

end
