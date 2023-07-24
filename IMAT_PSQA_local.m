function [G Gav]=IMAT_PSQA_local(DD,DTA,LDT)

d=dir();
%%pat=sprintf(d.name);
DD=DD/100;
LDT=LDT/100;

G1=zeros(length(d),1);
G2=zeros(length(d),1);
G_2p2mm5p_A1=zeros(length(d),1);
G_2p2mm5p_A2=zeros(length(d),1);

for i = 3:length(d)
    cd(d(i).name)
   
    startRow=263;
    endRow=303;
    d1=dir('*.txt');
    d2=dir('*.snc');

if length(d1)==2 && length(d2)==2   %%skip the file that doesn't contain txt and snc file

%Arc 1    
fileName1 = sprintf(d1(1).name);

  meas1 = importfile01(fileName1,startRow,endRow);

%Arc 2  
fileName2 = sprintf(d1(2).name);

  meas2 = importfile01(fileName2,startRow,endRow);  


    startRow=7;
    endRow=inf;

    %TPS 1
  fileName3 = sprintf(d2(1).name);
  tps1 = importfile02(fileName3,startRow,endRow);
 
  %TPs 2
   fileName4 = sprintf(d2(2).name);
   tps2 = importfile02(fileName4,startRow,endRow);
  
  %Normalization 
    meas1=meas1./max(meas1(:));
  tps1=tps1./max(tps1(:));
   meas2=meas2./max(meas2(:));
  tps2=tps2./max(tps2(:)); 
  
  g1=gamma_tool(meas1,tps1,DD,DTA,LDT,'local');
  G1(i)=g1;
  g2=gamma_tool(meas2,tps2,DD,DTA,LDT,'local');
  G2(i)=g2;
  ggA1=gamma_tool(meas1,tps1,0.02,2,0.05,'global');  %%For filtering purpose
  G_2p2mm5p_A1(i)=ggA1;
  ggA2=gamma_tool(meas2,tps2,0.02,2,0.05,'global');
  G_2p2mm5p_A2(i)=ggA2;
  
end

cd ..
  
end

G_arc1=G_2p2mm5p_A1;
G_arc2=G_2p2mm5p_A2;

f1=find(G_arc1>=0.935); %%Not all files are matched,  so we are just focussing on the files we want.
f2=find(G_arc2>=0.935);

G1=G1(f1)*100;
G2=G2(f2)*100;

G1av=mean(G1);
G2av=mean(G2);

G1_std=std(G1);
G2_std=std(G2);

G=[G1;G2];

G_std=std(G(:));
Gav=mean(G(:));
Gav=[Gav G_std];

DD=DD*100;
DTA=DTA;
LDT=LDT*100;

formatSpec1='G_local_%dp_%dmm_%dp_IMAT.txt';
formatSpec2='Gav_local_%dp_%dmm_%dp_IMAT.txt';

str1 = sprintf(formatSpec1,DD,DTA,LDT);
str2 = sprintf(formatSpec2,DD,DTA,LDT);

cd ..

cd results
cd IMAT_g_local

T1 = table(G);
writetable(T1,str1,'Delimiter',' ');

T2 = table(Gav);
writetable(T2,str2,'Delimiter',' ');


cd ..
cd ..

cd IMAT
end