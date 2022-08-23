function [Ms, Mw, EN, ENr, SYe, SXe, SE, SG, OPP, yr, xc] = Extrema_BSA(img)
% Made by Tuan Nguyen, Tsvetkov, 2019.04.22
% BLOCK-SEGMENT SEARCH OF LOCAL EXTREME IMAGES BASED ON ANALYSIS OF BRIGHTNESSES OF RELATED PIXELS AND AREAS
% BSA (Block Segment Analysis)
%---------------------------------
% Ms: segmented image of local extrema, [0, |SN|]
% Mw: segmented image for SRG
% SE: stack of number of local extrema
% SG: stack of extreme values
% SYe, SXe: stacks of coordinates of the first extreme pixels for Seeded Region Growing
% SN: label of processed regions
% img: input graysacle or RGB image
% BH: block matrix of maxima
% BL: block matrix of minima
% EC: Stack of number elements of segment
% ES: Stack of labels
% EN: number of extrema
% ENr: number of border pixels of extrema
% ENm: number of Maxima
% EM: label matrix of extrem regions
% ER: check matrix for RG
% SYr, SXr: coordinates, seed value of border pixels of extrema.
% OPP: Operations per Pixel
% Memory resource = (64)*n (n-the number of pixels of image of 8bit/pixel).
%
% Convert RGB to Grayscale image
if ndims(img)==3
    img = rgb2gray(img);
end
img = double(img);
%
img = padarray(img, [1 1],'replicate','both');%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input paramaters
[Y, X] = size(img);
Len = Y*X;
BH = ones(Y,X);  % Local Maxima Matrix =1
BL = ones(Y,X);  % Local Minima Matrix =1
Ms = zeros(Y,X); % Label matrix of extrema =8
Mw = zeros(Y,X); % Label matrix for RG
EM = zeros(Y,X); % Label matrix for RG =8
ER = zeros(Y,X); % Check matrix for RG =8
N=0; %number of operations
EN=0; %number of extrema
SE=zeros(Len,2); %Stack of counter of extrema and extreme value 8, 8
ENr=0; SYr=zeros(1,Len); SXr=zeros(1,Len); %Stacks of region pixels 8, 8
SYe=zeros(Len,1); SXe=zeros(Len,1);
SN=0;  SY=zeros(1,Len);  SX=zeros(1,Len); %Stacks for region growing -> 8,8(YX)
% Create grids of local extrema for approximaion
yr=zeros(1,Y-2); yr(1)=1; yr(end)=1; 
xc=zeros(1,X-2); xc(1)=1; xc(end)=1;
%
y=2; 
while y<Y
    x=2;
    while x<X        
        % Scanning Mask: [6 7 8; 2 0 1; 5 4 3]
        S6=img(y-1,x-1); S7=img(y-1,x); S8=img(y-1,x+1);
        S2=img(y,x-1);   S0=img(y,x);   S1=img(y,x+1);
        S5=img(y+1,x-1); S4=img(y+1,x); S3=img(y+1,x+1);
        N=N+11;
        % 1.Find Local Mixima        
        if BH(y,x)==1
            %N=N+1;
            if S0>=S1
                N=N+2;
                BH(y,x+1)=0;%BH(S1)
                if S0>=S2%
                    N=N+2;
                    BH(y,x-1)=0;%BH(S2)
                    if S0>=S3
                        N=N+2;
                        BH(y+1,x+1)=0;%BH(S3) 
                        if S0>=S4
                            N=N+2;
                            BH(y+1,x)=0;%BH(S4)
                            if S0>=S5
                                N=N+2;
                                BH(y+1,x-1)=0;%BH(S5)
                                if S0>=S6%
                                    N=N+2;
                                    BH(y-1,x-1)=0;%BH(S6)
                                    if S0>=S7%
                                        N=N+2;
                                        BH(y-1,x)=0;%BH(S7)
                                        if S0>=S8
                                            N=N+2;
                                            BH(y-1,x+1)=0;%BH(S8)
                                            %Find group local Maxima
                                            if (S0==S1)||(S0==S3)||(S0==S4)||(S0==S5)||(S0==S2)||(S0==S6)||(S0==S7)||(S0==S8) 
                                                SN=SN+1;
                                                coun=1; %Ceny=y; Cenx=x;
                                                EM(y,x)=SN; ER(y,x)=SN;
                                                SP=1; SY(SP)=y; SX(SP)=x;
                                                ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x; 
                                                BH(y,x)=0; BL(y,x)=0;
                                                Emax=0; Emin=0; N=N+14;
                                                while(SP>0)                                                    
                                                    y1=SY(SP); 
                                                    x1=SX(SP); 
                                                    SP=SP-1;
                                                    N=N+4;
                                                    %
                                                    i=-1;
                                                        while i<2
                                                            j=-1;
                                                            while j<2
                                                                if (j==0)&&(i==0)
                                                                    N=N+1;
                                                                else
                                                                    N=N+2;%!
                                                                    ny=y1+i;
                                                                    nx=x1+j;
                                                                    if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&ER(ny,nx)~=SN
                                                                        ER(ny,nx)=SN; N=N+1;
                                                                        if img(ny,nx)==img(y1,x1)&&EM(ny,nx)==0                              
                                                                            EM(ny,nx)=SN; %Ceny=Ceny+ny; Cenx=Cenx+nx;
                                                                            BH(ny,nx)=0; BL(ny,nx)=0;
                                                                            SP=SP+1; SY(SP)=ny; SX(SP)=nx;
                                                                            if (ny==1) || (ny==Y) || (nx==1) || (nx==X)
                                                                                N=N+1;
                                                                            else
                                                                                ENr=ENr+1; SYr(ENr)=ny; SXr(ENr)=nx;
                                                                                coun=coun+1;
                                                                                N=N+4;
                                                                            end
                                                                            N=N+6;
                                                                        else
                                                                            N=N+1;
                                                                            if img(ny,nx)>img(y1,x1) %min
                                                                                Emin=1; 
                                                                                BL(ny,nx)=0;
                                                                                N=N+3;
                                                                            else %max
                                                                                Emax=1;
                                                                                BH(ny,nx)=0; 
                                                                                N=N+2;
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                                j=j+1; 
                                                            end
                                                            i=i+1; 
                                                        end
                                                        %
                                                end
                                                %
                                                if (Emax==1)&&(Emin==1) %non-extremum
                                                    ENr=ENr-coun;
                                                    N=N+2;
                                                else
                                                    EN=EN+1;
                                                    SYe(EN)=y;%SYr(ENr-floor(coun/1.5));%floor(Ceny/coun); 
                                                    SXe(EN)=x;%SXr(ENr-floor(coun/1.5));%floor(Cenx/coun);
                                                    yr(y-1)=1; xc(x-1)=1;
%                                                     if (y>3)&&(y<Y-3)&&(yr(y-2)+yr(y)==0)
%                                                         yr(y-1)=1;
%                                                     end
%                                                     if (x>3)&&(x<X-3)&&(xc(x-2)+xc(x)==0)
%                                                         xc(x-1)=1;
%                                                     end
                                                    NS=EN*(Emax-Emin);
                                                    Mw(y,x)=NS;
                                                    SE(EN,1)=coun;
                                                    SE(EN,2)=round(img(y,x));
                                                    for j=(ENr-coun+1):ENr
                                                        y1=SYr(j); x1=SXr(j);
                                                        Ms(y1,x1)=NS;
                                                    end
                                                    N=N+7+coun*3;%!
                                                    %
                                                end
                                                %
                                            else % Find single local Maxima
                                                EN=EN+1;
                                                SYe(EN)=y; SXe(EN)=x; Mw(y,x)=EN;
                                                Ms(y,x)=EN; BH(y,x)=0; BL(y,x)=0;
                                                SE(EN,1)=1;
                                                SE(EN,2)=round(img(y,x));
                                                ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x; yr(y-1)=1; xc(x-1)=1;
                                                N=N+12;
                                            end
                                             %
                                        else %S0<S8
                                            N=N+2;
                                            BL(y-1,x+1)=0;
                                            BH(y,x)=0;%BH(S0)
                                        end
                                    else %S0<S7
                                        N=N+2;
                                        BL(y-1,x)=0;
                                        BH(y,x)=0;%BH(S0)
                                    end
                                else %S0<S6
                                    N=N+2;
                                    BL(y-1,x-1)=0;
                                    BH(y,x)=0;%BH(S0)
                                end
                            else %S0<S5
                                N=N+2;
                                BL(y+1,x-1)=0;
                                BH(y,x)=0;
                            end
                        else %S0<S4
                            N=N+2;
                            BL(y+1,x)=0;%BL(S4)
                            BH(y,x)=0;%BH(S0)
                        end
                    else %S0<S3
                        N=N+2;
                        BL(y+1,x+1)=0;%BL(S3)
                        BH(y,x)=0;%BH(S0)
                    end
                else %S0<S2
                    N=N+2;
                    BL(y,x-1)=0;%BL(S2)
                    BH(y,x)=0;%BH(S0)
                end
            else %S0<S1
                N=N+2;
                BL(y,x+1)=0;%BL(S1)
                BH(y,x)=0;%BH(S0)
            end
        end
       % 2.Find Local Mimima
        if BL(y,x)==1
            %N=N+1;
            if S0<=S1
                N=N+2;
                BL(y,x+1)=0;%BL(S1)
                if S0<=S2
                    N=N+2;
                    BL(y,x-1)=0;%BL(S2)
                    if S0<=S3
                        N=N+2;
                        BL(y+1,x+1)=0;%BL(S3)
                        if S0<=S4
                            N=N+2;
                            BL(y+1,x)=0;%BL(S4)
                            if S0<=S5
                                N=N+2;
                                BL(y+1,x-1)=0;%BL(S5)
                                if S0<=S6
                                    N=N+2;
                                    BL(y-1,x-1)=0;%BL(S6)
                                    if S0<=S7
                                        N=N+2;
                                        BL(y-1,x)=0;%BL(S7)
                                        if S0<=S8
                                            N=N+2;
                                            BL(y-1,x+1)=0;%BL(S8)
                                            %Find group local Minima
                                            if (S0==S1)||(S0==S3)||(S0==S4)||(S0==S5)||(S0==S2)||(S0==S6)||(S0==S7)||(S0==S8)     
                                                SN=SN+1;
                                                coun=1; %Ceny=y; Cenx=x; 
                                                EM(y,x)=SN; ER(y,x)=SN;
                                                SP=1; SY(SP)=y; SX(SP)=x;
                                                ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x;
                                                BH(y,x)=0; BL(y,x)=0;%!
                                                Emax=0; Emin=0; N=N+14;
                                                while(SP>0)                                                    
                                                    y1=SY(SP); 
                                                    x1=SX(SP); 
                                                    SP=SP-1; 
                                                    N=N+4;
                                                    %
                                                    i=-1;
                                                        while i<2
                                                            j=-1;
                                                            while j<2
                                                                if (j==0)&&(i==0)
                                                                    N=N+1;
                                                                else
                                                                    N=N+2;%!
                                                                    ny=y1+i;
                                                                    nx=x1+j;
                                                                    if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&ER(ny,nx)~=SN
                                                                        ER(ny,nx)=SN; N=N+1;%!
                                                                        if img(ny,nx)==img(y1,x1)&&EM(ny,nx)==0                              
                                                                            EM(ny,nx)=SN; %Ceny=Ceny+ny; Cenx=Cenx+nx;
                                                                            BH(ny,nx)=0; BL(ny,nx)=0;%!
                                                                            SP=SP+1; SY(SP)=ny; SX(SP)=nx;
                                                                            if (ny==1) || (ny==Y) || (nx==1) || (nx==X)
                                                                                N=N+1;
                                                                            else
                                                                                ENr=ENr+1; SYr(ENr)=ny; SXr(ENr)=nx;
                                                                                coun=coun+1;
                                                                                N=N+4;
                                                                            end
                                                                            N=N+6;
                                                                        else
                                                                            N=N+1;
                                                                            if img(ny,nx)<img(y1,x1) %max
                                                                                Emax=1; 
                                                                                BH(ny,nx)=0;
                                                                                N=N+3;
                                                                            else %min
                                                                                Emin=1;
                                                                                BL(ny,nx)=0; 
                                                                                N=N+2;
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                                j=j+1; 
                                                            end
                                                            i=i+1;
                                                        end
                                                        %
                                                end
                                                % Flag
                                                if (Emax==1)&&(Emin==1) %non-extremum
                                                    ENr=ENr-coun;
                                                    N=N+2;
                                                else
                                                    EN=EN+1; 
                                                    SYe(EN)=y;%SYr(ENr-floor(coun/1.5));%floor(Ceny/coun); 
                                                    SXe(EN)=x;%SXr(ENr-floor(coun/1.5));%floor(Cenx/coun); 
                                                    yr(y-1)=1; xc(x-1)=1;
%                                                     if (y>3)&&(y<Y-3)&&(yr(y-2)+yr(y)==0)
%                                                         yr(y-1)=1;
%                                                     end
%                                                     if (x>3)&&(x<X-3)&&(xc(x-2)+xc(x)==0)
%                                                         xc(x-1)=1;
%                                                     end
                                                    NS=EN*(Emax-Emin);
                                                    Mw(y,x)=NS;
                                                    SE(EN,1)=coun;
                                                    SE(EN,2)=round(img(y,x));
                                                    for j=(ENr-coun+1):ENr
                                                        y1=SYr(j); x1=SXr(j);
                                                        Ms(y1,x1)=NS;%!
                                                    end
                                                    N=N+7+coun*3;%!
                                                    %
                                                end
                                                %
                                            else %Find single local Minima
                                                EN=EN+1; SYe(EN)=y; SXe(EN)=x; Mw(y,x)=-EN;
                                                Ms(y,x)=-EN; BL(y,x)=0; BH(y,x)=0;
                                                SE(EN,1)=1;
                                                SE(EN,2)=round(img(y,x));
                                                ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x; yr(y-1)=1; xc(x-1)=1;
                                                N=N+12;
                                            end
                                            %
                                        else %S0>S8
                                            N=N+2;
                                            BH(y-1,x+1)=0;%BH(S8)
                                            BL(y,x)=0;%BL(S0)
                                        end
                                    else %S0>S7
                                        N=N+2;
                                        BH(y-1,x)=0;%BH(S7)
                                        BL(y,x)=0;%BL(S0)
                                    end
                                else %S0>S6
                                    N=N+2;
                                    BH(y-1,x-1)=0;%BH(S6)
                                    BL(y,x)=0;%BL(S0)
                                end
                            else %S0>S5
                                N=N+2;
                                BH(y+1,x-1)=0;%BH(S5)
                                BL(y,x)=0;%BL(S0)
                            end
                        else %S0>S4
                            N=N+2;
                            BH(y+1,x)=0;%BH(S4)
                            BL(y,x)=0;%BL(S0)
                        end
                    else %S0>S3
                        N=N+2;
                        BH(y+1,x+1)=0;%BH(S3)
                        BL(y,x)=0;%BL(S0)
                    end
                else %S0>=S2
                    N=N+2;
                    BH(y,x-1)=0;%BH(S2)
                    BL(y,x)=0;%BL(S0)
                end
            else %S0>S1
                N=N+2;
                BH(y,x+1)=0;%BH(S1)
                BL(y,x)=0;%BL(S0)
            end
        end     
        x = x+1;
    end
    y = y+1;
end
%
OPP = N/(Y*X);
%
SE(EN+1:end,:)=[];
SG=SE; SG(:,1)=1;
%
SYe(1:EN)=SYe(1:EN)-1; SYe((y-2)*(X-2)+1:end)=[];
SXe(1:EN)=SXe(1:EN)-1; SXe((y-2)*(X-2)+1:end)=[];
%SYr(1:ENr)=SYr(1:ENr)-1; SXr(1:ENr)=SXr(1:ENr)-1;
%
Ms = Ms(2:Y-1,2:X-1);
Mw = Mw(2:Y-1,2:X-1);
%
Y=Y-2; X=X-2;
%
end