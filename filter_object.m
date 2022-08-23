function [imF, imS, SMI2, NSS, SYs, SXs] = filter_object(img, SMI, SYe, SXe, SA, SE, EN, SIZE, T1, T2)
% Object filtering by size, made by Tuan Nguyen, 28.05.2021
if ndims(img)==3
    img = rgb2gray(img);
end
img = double(img);
[Y,X]=size(img);
SMI2=zeros(Y,X);
imS=zeros(Y,X);
imF=img;
SY=zeros(1,Y*X);  SX=zeros(1,Y*X);
NSS=0; SYs=zeros(1,Y*X); SXs=zeros(1,Y*X);
%1. Object localization
for k=1:EN
    if  (SE(k,2)>=T1)&&(SE(k,2)<=T2) || (SA(k,1)<=SIZE)
        y=SYe(k); x=SXe(k);
        SN=SMI(y,x);
        NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
        imS(y,x)=SN; SMI2(y,x)=1;
        SP=1; SY(SP)=y; SX(SP)=x;
        while SP>0
            y1=SY(SP); 
            x1=SX(SP); 
            SP=SP-1;
            for j=-1:1
                for i=-1:1
                    if (j==0)&&(i==0)
                    else
                        ny=y1+j; nx=x1+i;
                        if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&imS(ny,nx)==0
                            if SMI(ny,nx)==SMI(y1,x1)
                                SMI2(ny,nx)=1;
                                imS(ny,nx)=SN;
                                SP=SP+1; SY(SP)=ny; SX(SP)=nx;
                                NSS=NSS+1; SYs(NSS)=ny; SXs(NSS)=nx;
                            end
                        end
                    end
                end
            end
        end
    end
end

%2. Small object Deleting
k=1;
while k<=NSS && NSS<=Y*X
    y=SYs(k); x=SXs(k); AVG=0; idx=0;
    for j=-1:1
        for i=-1:1
            if (j==0)&&(i==0)
            else
                ny=y+j; nx=x+i;
                if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&SMI2(ny,nx)==0
                    AVG=AVG+imF(ny,nx); idx=idx+1;
                end
            end
        end
    end
    if idx==0
        NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
    else
        imF(y,x)=round(AVG/idx); SMI2(y,x)=0;
    end
    k=k+1;
end
%
end
%%
% function [imF, imS, SMI2, NSS, SYs, SXs] = filter_object(img, SMI, SYe, SXe, SA, SE, EN, SIZE1, SIZE2, T1, T2)
% % Object filtering by size, made by Tuan Nguyen, 28.05.2021
% if ndims(img)==3
%     img = rgb2gray(img);
% end
% img = double(img);
% [Y,X]=size(img);
% SMI2=zeros(Y,X);
% imS=zeros(Y,X);
% imF=img;
% SY=zeros(1,Y*X);  SX=zeros(1,Y*X);
% NSS=0; SYs=zeros(1,Y*X); SXs=zeros(1,Y*X);
% %1. Object localization
% for k=1:EN
%     if  (SE(k,2)>=T1)&&(SE(k,2)<=T2)&&(SA(k,1)<=SIZE2) %|| (SE(k,2)==SA(k,2)) %%|| (SE(k,2)>=T2) %(SA(k,1)<=SIZE1)  || (SE(k,2)>=T1)&&(SE(k,2)<=T2)&&(SA(k,1)<=SIZE2) ||
%         y=SYe(k); x=SXe(k);
%         SN=SMI(y,x);
%         NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
%         imS(y,x)=SN; SMI2(y,x)=1;
%         SP=1; SY(SP)=y; SX(SP)=x;
%         while SP>0
%             y1=SY(SP); 
%             x1=SX(SP); 
%             SP=SP-1;
%             for j=-1:1
%                 for i=-1:1
%                     if (j==0)&&(i==0)
%                     else
%                         ny=y1+j; nx=x1+i;
%                         if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&imS(ny,nx)==0&&SMI(ny,nx)==SMI(y1,x1)
%                             SMI2(ny,nx)=1;
%                             imS(ny,nx)=SN;
%                             SP=SP+1; SY(SP)=ny; SX(SP)=nx;
%                             NSS=NSS+1; SYs(NSS)=ny; SXs(NSS)=nx;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% %2. Small object Deleting
% k=1;
% while k<=NSS && NSS<=Y*X
%     y=SYs(k); x=SXs(k); AVG=0; idx=0;
%     for j=-1:1
%         for i=-1:1
%             if (j==0)&&(i==0)
%             else
%                 ny=y+j; nx=x+i;
%                 if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&SMI2(ny,nx)==0%&&(imF(ny,nx)<T1)&&(imF(ny,nx)>T2)
%                     AVG=AVG+imF(ny,nx); idx=idx+1;
%                 end
%             end
%         end
%     end
%     if idx==0
%         NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
%     else
%         imF(y,x)=round(AVG/idx); SMI2(y,x)=0;
%     end
%     k=k+1;
% end
% % for k=1:NSS
% %     y=SYs(k); x=SXs(k); label=SMI(y,x); idx=abs(label); AVG=round(SA(idx,2));
% %     imF(y,x)=AVG;
% % end
% %
% end
%%
% function [imF, imS, SMI2, SM, ENr, SYr, SXr, NSS, SA] = filter_object(img, SMI, SYe, SXe, SA, EN, SIZE)
% % Object filtering by size, made by Tuan Nguyen, 28.05.2021
% if ndims(img)==3
%     img = rgb2gray(img);
% end
% img = double(img);
% [Y,X]=size(img);
% SMI2=zeros(Y,X);
% imS=zeros(Y,X);
% imF=img;
% SY=zeros(1,Y*X);  SX=zeros(1,Y*X);
% ENr=0; SYr=zeros(1,Y*X); SXr=zeros(1,Y*X);
% %1. Object localization
% for k=1:EN
%     if (SA(k,1)>0)&&(SA(k,1)<=SIZE)
%         y=SYe(k); x=SXe(k); label=SMI(y,x); index=abs(label); SA(index,1)=0; SA(index,2)=0;
%         SN=SMI(y,x);
%         imS(y,x)=SN; SMI2(y,x)=1;
%         SP=1; SY(SP)=y; SX(SP)=x;
%         while SP>0
%             y1=SY(SP); 
%             x1=SX(SP); 
%             SP=SP-1;
%             for j=-1:1
%                 for i=-1:1
%                     if (j==0)&&(i==0)
%                     else
%                         ny=y1+j; nx=x1+i;
%                         if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&imS(ny,nx)==0&&SMI(ny,nx)==SMI(y1,x1)
%                             SMI2(ny,nx)=1;
%                             imS(ny,nx)=SN;
%                             SP=SP+1; SY(SP)=ny; SX(SP)=nx;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% %2.Object region merging
% SM=zeros(Y,X);
% SYs=zeros(1,Y*X); 
% SXs=zeros(1,Y*X);
% NS=0; NSS=0;
% y=0;
% while y<Y
%     x=0;
%     while x<X
%         if (SM(y+1,x+1)==0)&&(SMI2(y+1,x+1)==1)
%             NS=NS+1; 
%             NSS=NSS+1; SYs(NSS)=y+1; SXs(NSS)=x+1;
%             SM(y+1,x+1)=NS; 
%             SP=1;
%             SY(SP)=y;
%             SX(SP)=x;
%             while (SP>0)
%                 y1=SY(SP);
%                 x1=SX(SP);
%                 SP=SP-1;
%                 for j=-1:1
%                     for i=-1:1
%                         if (j==0)&&(i==0)
%                         else
%                             if (y1+j+1>0)&&(x1+i+1>0)&&(y1+j+1<=Y)&&(x1+i+1<=X)&&(SMI2(y1+j+1,x1+i+1)==1)&&(SM(y1+j+1,x1+i+1)==0)
%                                 SM(y1+j+1,x1+i+1)=NS;
%                                 SP=SP+1; SY(SP)=y1+j; SX(SP)=x1+i;
%                                 NSS=NSS+1; SYs(NSS)=y1+j+1; SXs(NSS)=x1+i+1;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         x=x+1;
%     end
%     y=y+1;
% end
% %3. Small object Deleting
% SMI3=SMI2;
% k=1;
% while k<=NSS
%     y=SYs(k); x=SXs(k); AVG=0; idx=0;
%     for j=-1:1
%         for i=-1:1
%             if (j==0)&&(i==0)
%             else
%                 ny=y+j; nx=x+i;
%                 if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&(SMI3(ny,nx)==0)
%                     AVG=AVG+img(ny,nx); idx=idx+1;
%                 end
%             end
%         end
%     end
%     if idx==0
%         NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
%     else
%         imF(y,x)=round(AVG/idx); SMI3(y,x)=0;
%     end
%     k=k+1;
% end
% %
% end
%%