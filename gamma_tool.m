function g = gamma_tool(Ref,Eval,DD,DTA,threshold,localorglobal)

%DTA=Distance-to-agreement, units depend on the dose distribution spatial domain
%DD=Dose Difference, specified in decimals

%threshold is the minimum percentage of low dose threshold, specified in decimals

%Obtain the sizes of both dose distributions
s_Ref = size(Ref);
s_Eval = size(Eval);


%Correct registration error between Arccheck data and TPS

MD=zeros(s_Eval);

if s_Ref~=s_Eval
for i=1:s_Ref(1)
    for j=1:s_Ref(2)
        if Ref(i,j)==0
            continue
        end
        
       
        MD(((i+1)*5)+1,(j+2)*5)=Ref(i,j);  %%Resize Ref dose and put the diode dose value into it
        
    end
end
Ref=MD;
end



%%Declare again
s_Ref = size(Ref);
s_Eval = size(Eval);

Gamma = zeros(s_Ref);

threshold = threshold*max(Ref( : )) ;

searchRange = round(3*DTA+1);

%Create the DTAsearch range and the DTA_Test
%The DTA_Test is the same for each evaluated pixel
DTA_Test = zeros(2*searchRange+1 ,2*searchRange+1) ;

%Sets the origin position to 1
DTA_Test(searchRange+1, searchRange+1)=1;

%Calculates the euclidean distance to all points from the origin and divides with DTA and quadrates this matrice.
DTA_Test =(bwdist(DTA_Test)/DTA).^2;
 
    Ref_temp = zeros(s_Ref+searchRange*2);
    Ref_temp(searchRange+1:s_Ref(1)+searchRange,searchRange+1:s_Ref(2)+searchRange) = Ref;
    Eval_temp = zeros(s_Eval+searchRange*2) ;
    Eval_temp(searchRange+1: s_Ref(1)+searchRange, searchRange+1:s_Ref(2)+searchRange) = Eval;
    Ref = Ref_temp ;
    Eval = Eval_temp ;
    Gamma = zeros(size(Ref));

    if strcmpi(localorglobal,'global')
        DD = zeros(size(Ref))+DD*max(Ref(:));
    elseif strcmpi (localorglobal,'local')
        DD = Ref*DD;
    end


x=-(searchRange):1:(searchRange);
y=-(searchRange):1:(searchRange);

if s_Ref == s_Eval
    for i = searchRange+1:s_Ref(1)+searchRange
        for j = searchRange+1:s_Ref(2)+searchRange
            
            if Ref(i,j)<threshold %To ignore unneccesary pixels
            Gamma(i,j)= 0 ;
            continue
            end

            DD_function = ((Ref(i,j)-Eval(i+x,j+y))) ;
            DD_Test = (DD_function./DD(i,j)).^2 ;
            %The Gamma matrice
            G = DTA_Test+ DD_Test;
            %minimum of Gamma
            Gamma(i,j)=min(G(:));
        end
    end



end

Gamma = sqrt(Gamma);

%gc=length(find(Gamma<=1))-length(find(Gamma==0));
g=(length(find(Gamma<=1))-length(find(Gamma==0)))/length(find(Gamma~=0));


if s_Ref==s_Eval
Gamma(1 : searchRange , :) = [ ];
Gamma(: , 1 : searchRange) = [ ];
Gamma(s_Ref(1)+1: s_Ref(1)+searchRange , :) = [ ];
Gamma(: , s_Ref(2)+1: s_Ref(2)+searchRange) = [ ];

end



%figure
%imagesc(Gamma)
%colormap ('jet')
%colorbar
%caxis([0 2])

end