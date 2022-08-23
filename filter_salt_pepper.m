function [imF, imS, SME2, NSS, SYs, SXs] = filter_salt_pepper(imnoise, SME, SYe, SXe, SE, EN, T1, T2)
% Salt&Pepper Filtering made by Tuan Nguyen, 28.05.2021
% using thresholds (T1=0) and (T2=2^B-1) for filtering.
% we can use two thresolds for salt and pepper: <=T1 and >=T2
if ndims(imnoise)==3
    imnoise = rgb2gray(imnoise);
end
imnoise = double(imnoise);
[Y,X]=size(imnoise);
SME2=zeros(Y,X);
imS=zeros(Y,X);
imF=imnoise;
SY=zeros(1,Y*X);  SX=zeros(1,Y*X);
NSS=0; SYs=zeros(1,Y*X); SXs=zeros(1,Y*X);
%1.Salt & Pepper localization
for k=1:EN
    if SE(k,2)<=T1 || SE(k,2)>=T2
        y=SYe(k); x=SXe(k);
        SN=SME(y,x);
        imS(y,x)=SN; 
        SME2(y,x)=1;
        SP=1; SY(SP)=y; SX(SP)=x;
        NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
        while SP>0
            y1=SY(SP); 
            x1=SX(SP); 
            SP=SP-1;
            for j=-1:1
                for i=-1:1
                    if (j==0)&&(i==0)
                    else
                        ny=y1+j; nx=x1+i;
                        if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&imS(ny,nx)==0&&SME(ny,nx)==SME(y1,x1)
                            SME2(ny,nx)=1;
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
%2. Salt&Pepper Deleting
k=1;
while k<=NSS
    y=SYs(k); x=SXs(k); AVG=0; idx=0;
    for j=-1:1
        for i=-1:1
            if (j==0)&&(i==0)
            else
                ny=y+j; nx=x+i;
                if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&SME2(ny,nx)==0%(imF(ny,nx)>T1)&&(imF(ny,nx)<T2)
                    AVG=AVG+imF(ny,nx); idx=idx+1;
                end
            end
        end
    end
    if idx==0
        NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
    else
        imF(y,x)=round(AVG/idx); SME2(y,x)=0;
    end
    k=k+1;
end
%
end
%%
%Удаление шума соли и перца с изображений
%http://wiki.technicalvision.ru/index.php/%D0%A4%D0%B8%D0%BB%D1%8C%D1%82%D1%80%D0%B0%D1%86%D0%B8%D1%8F_%D0%B1%D0%B8%D0%BD%D0%B0%D1%80%D0%BD%D1%8B%D1%85_%D0%B8%D0%B7%D0%BE%D0%B1%D1%80%D0%B0%D0%B6%D0%B5%D0%BD%D0%B8%D0%B9
%https://russianblogs.com/article/5382204114/
%http://jre.cplire.ru/jre/apr09/7/text.html
%https://medium.com/analytics-vidhya/remove-salt-and-pepper-noise-with-median-filtering-b739614fe9db
%https://en.wikipedia.org/wiki/Salt-and-pepper_noise
%https://www.mathworks.com/help/vision/ug/remove-salt-and-pepper-noise-from-images.html
%https://coderoad.ru/12452942/%D0%A3%D0%B4%D0%B0%D0%BB%D0%B5%D0%BD%D0%B8%D0%B5-%D0%BF%D0%B5%D1%80%D1%86%D0%BE%D0%B2%D0%BE%D0%B3%D0%BE-%D1%88%D1%83%D0%BC%D0%B0-%D0%B8-%D1%81%D0%BE%D0%BB%D0%B5%D0%BD%D0%BE%D0%B3%D0%BE-%D1%88%D1%83%D0%BC%D0%B0-%D0%B8%D0%B7-%D0%B8%D0%B7%D0%BE%D0%B1%D1%80%D0%B0%D0%B6%D0%B5%D0%BD%D0%B8%D0%B9
%https://www.google.com/search?q=how+to+Remove+Salt+and+Pepper+Noise+from+Image&oq=how+to+Remove+Salt+and+Pepper+Noise+from+Image&aqs=chrome..69i57j0i19i22i30l3.7930j0j1&sourceid=chrome&ie=UTF-8
%% 01
% function [imF, imS, SME2, SM, ENr, SYr, SXr, NS, SE, NSS] = filter_salt_pepper(img, SME, SYe, SXe, SE, EN, T1, T2)
% % Salt&Pepper Filtering made by Tuan Nguyen, 28.05.2021
% % using thresholds (T1=0) and (T2=255) for filtering.
% % we can use two thresolds for salt and pepper: <=T1 and >=T2
% if ndims(img)==3
%     img = rgb2gray(img);
% end
% img = double(img);
% [Y,X]=size(img);
% SME2=zeros(Y,X);
% imS=zeros(Y,X);
% imF=img;
% SY=zeros(1,Y*X);  SX=zeros(1,Y*X);
% ENr=0; SYr=zeros(1,Y*X); SXr=zeros(1,Y*X);
% %1.Salt & Pepper localization
% for k=1:EN
%     if SE(k,2)<=T1 || SE(k,2)>=T2
%         y=SYe(k); x=SXe(k);
%         SN=SME(y,x);
%         imS(y,x)=SN; 
%         SME2(y,x)=1;
%         SP=1; SY(SP)=y; SX(SP)=x;
%         ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x;
%         while SP>0
%             y1=SY(SP); 
%             x1=SX(SP); 
%             SP=SP-1;
%             for j=-1:1
%                 for i=-1:1
%                     if (j==0)&&(i==0)
%                     else
%                         ny=y1+j; nx=x1+i;
%                         if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&imS(ny,nx)==0&&SME(ny,nx)==SME(y1,x1)
%                             SME2(ny,nx)=1;
%                             imS(ny,nx)=SN;
%                             SP=SP+1; SY(SP)=ny; SX(SP)=nx;
%                             ENr=ENr+1; SYr(ENr)=ny; SXr(ENr)=nx;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% %2.Salt&Pepper region merging
% SM=zeros(Y,X);
% SYs=zeros(1,Y*X); 
% SXs=zeros(1,Y*X);
% NS=0; NSS=0;
% for k=1:ENr
%     y=SYr(k); x=SXr(k);
%     if SM(y,x)==0
%         NS=NS+1; 
%         NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
%         SM(y,x)=NS; 
%         SP=1;
%         SY(SP)=y;
%         SX(SP)=x;
%         while (SP>0)
%             y1=SY(SP);
%             x1=SX(SP);
%             SP=SP-1;
%             for j=-1:1
%                 for i=-1:1
%                     if (j==0)&&(i==0)
%                     else
%                         ny=y1+j; nx=x1+i;
%                         if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&(SME2(ny,nx)==1)&&(SM(ny,nx)==0)
%                             SM(ny,nx)=NS;
%                             SP=SP+1; SY(SP)=ny; SX(SP)=nx;
%                             NSS=NSS+1; SYs(NSS)=ny; SXs(NSS)=nx;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% %3. Salt&Pepper Deleting
% k=1;
% while k<=NSS
%     y=SYs(k); x=SXs(k); AVG=0; idx=0;
%     for j=-1:1
%         for i=-1:1
%             if (j==0)&&(i==0)
%             else
%                 ny=y+j; nx=x+i;
%                 if (ny>0)&&(nx>0)&&(ny<=Y)&&(nx<=X)&&(imF(ny,nx)>T1)&&(imF(ny,nx)<T2)
%                     AVG=AVG+imF(ny,nx); idx=idx+1;
%                 end
%             end
%         end
%     end
%     if idx==0
%         NSS=NSS+1; SYs(NSS)=y; SXs(NSS)=x;
%     else
%         imF(y,x)=round(AVG/idx);
%     end
%     k=k+1;
% end
% %
% end
%% Удаление шума соли и перца с изображений
%http://wiki.technicalvision.ru/index.php/%D0%A4%D0%B8%D0%BB%D1%8C%D1%82%D1%80%D0%B0%D1%86%D0%B8%D1%8F_%D0%B1%D0%B8%D0%BD%D0%B0%D1%80%D0%BD%D1%8B%D1%85_%D0%B8%D0%B7%D0%BE%D0%B1%D1%80%D0%B0%D0%B6%D0%B5%D0%BD%D0%B8%D0%B9
%https://russianblogs.com/article/5382204114/
%http://jre.cplire.ru/jre/apr09/7/text.html
%https://medium.com/analytics-vidhya/remove-salt-and-pepper-noise-with-median-filtering-b739614fe9db
%https://en.wikipedia.org/wiki/Salt-and-pepper_noise
%https://www.mathworks.com/help/vision/ug/remove-salt-and-pepper-noise-from-images.html
%https://coderoad.ru/12452942/%D0%A3%D0%B4%D0%B0%D0%BB%D0%B5%D0%BD%D0%B8%D0%B5-%D0%BF%D0%B5%D1%80%D1%86%D0%BE%D0%B2%D0%BE%D0%B3%D0%BE-%D1%88%D1%83%D0%BC%D0%B0-%D0%B8-%D1%81%D0%BE%D0%BB%D0%B5%D0%BD%D0%BE%D0%B3%D0%BE-%D1%88%D1%83%D0%BC%D0%B0-%D0%B8%D0%B7-%D0%B8%D0%B7%D0%BE%D0%B1%D1%80%D0%B0%D0%B6%D0%B5%D0%BD%D0%B8%D0%B9
%https://www.google.com/search?q=how+to+Remove+Salt+and+Pepper+Noise+from+Image&oq=how+to+Remove+Salt+and+Pepper+Noise+from+Image&aqs=chrome..69i57j0i19i22i30l3.7930j0j1&sourceid=chrome&ie=UTF-8
%%