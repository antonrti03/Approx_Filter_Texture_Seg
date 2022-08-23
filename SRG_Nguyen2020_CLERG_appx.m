%%
function [Ms, SP, SA, LOCAL_THRESHOLD_MIN, EN2, NODE, NGR] = SRG_Nguyen2020_CLERG_full(Ms, SYe, SXe, SG, EN, img, a, b, c, d)
% Segmentation algorithm CLERG+LE (proposed by Tuan Nguyen, Tsviatkou, 2020)
% Name: CLERG+LE full (Model CLERG+LE)
% Input:
% - Ms: segmented iamge;
% - SG: stack of one pixel extreman and their value 
% - SYe, SXe: stacks of coordinates of the first extreme pixels for Seeded Region Growing
% - SE: intensity of local extrema
% - EN: number of extrema
% - img: input grayscale image
% Output:
% - Ms: segmented image of one pixel extrema.
% - SP: number of loops for speed stack.
% - LOCAL_THRESHOLD_MIN: threshold
% - SPEED: speed stack.
%
% Convert RGB to grayscale image
if ndims(img)==3
    img = rgb2gray(img);
end
[Y, X] = size(img);
img = round(double(img));
%
NODE = sign(abs(Ms)); SXE=SXe; SYE=SYe;
NODE(1,1)=1; NODE(end:end)=1; NODE(1,end)=1; NODE(end,1)=1;
NGR  = img.*NODE;
% NGR(1,1)=img(1,1); 
% NGR(1,end)=img(1,end); 
% NGR(end,1)=img(end,1); 
% NGR(end,end)=img(end,end);
%
SA=SG;
%
maxx=max(max(double(img)));%B=8; 
GLOBAL_THRESHOLD = maxx; 
LOCAL_THRESHOLD_MIN=0;
loop=1; 
k=1;
EN2=1; SP=0;
while EN2>0
    EN2=0; LOCAL_THRESHOLD=maxx;
    while k<=EN
            y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label); EXT=SG(index,2); sgn=sign(label);
            id1=0; id2=0;
            for j=-1:1
                for i=-1:1
                    if (j==0)&&(i==0)
                    else
                        y1=y+j; x1=x+i;
                        if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
                            id1=id1+1; %avg=SA(index,2);
                            if img(y1,x1)==EXT-LOCAL_THRESHOLD_MIN*sgn
                                Ms(y1,x1)=label; 
                                EN=EN+1; SYe(EN)=y1; SXe(EN)=x1;
                                idx=SA(index,1);
                                %SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
                                SA(index,1) = idx+1;
                                XE=SXE(index); YE=SYE(index);
                                R = round(sqrt((x1-XE)^2+(y1-YE)^2));
                                %(12-24,17-33)
                                %(10-20,15-29)
                                %(8-16,13-25)
                                %(6-12,9-17)
                                %(4-8,7-13)
                                if mod(abs(y1-YE),a)==0&&mod(abs(x1-XE),a)==0&&R<=b || mod(abs(y1-YE),c)==0&&mod(abs(x1-XE),c)==0&&R>d
                                    NODE(y1,x1)=1; NGR(y1,x1)=img(y1,x1);
                                end
                                id2=id2+1;
                            else
                                if abs(EXT-img(y1,x1))>LOCAL_THRESHOLD_MIN && abs(EXT-img(y1,x1))<LOCAL_THRESHOLD
                                    LOCAL_THRESHOLD=abs(EXT-img(y1,x1));
                                end   
                            end
                        end
                    end 
                end
            end
            %
        if id1>id2
            EN2=EN2+1;
            SYe(EN2)=y; SXe(EN2)=x;
        end
        %
        k=k+1;
    end
    %
    if (EN2==0) || (LOCAL_THRESHOLD_MIN>=GLOBAL_THRESHOLD) || (LOCAL_THRESHOLD==maxx)
        break;
    else
        k=1; EN=EN2;
        LOCAL_THRESHOLD_MIN=LOCAL_THRESHOLD;
        loop=loop+1;
    end
end
%
% For Full Segmentation using Adams's SRG, 1994 for the last pixels
if EN2~=0
    EN=EN2;
    LOCAL_THRESHOLD_MIN=0; 
    k=1;
    EN2=1;
    while EN2>0
        EN2=0; LOCAL_THRESHOLD=maxx;
        while k<=EN
                y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label);
                id1=0; id2=0; 
                for j=-1:1
                    for i=-1:1
                        if (j==0)&&(i==0)
                        else
                            y1=y+j; x1=x+i;
                            if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
                                id1=id1+1; avg=SA(index,2);
                                if abs(img(y1,x1)-avg) <= LOCAL_THRESHOLD_MIN    
                                    Ms(y1,x1)=label; 
                                    EN=EN+1; SYe(EN)=y1; SXe(EN)=x1;
                                    idx=SA(index,1); 
                                    SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
                                    SA(index,1) = idx+1;
                                    id2=id2+1;
                                    %(10-20,15-29)
                                    %(8-16,13-25)
                                    %(6-12,9-17)
                                    %(4-8,7-13)
                                    %if mod(abs(y1-YE),4)==0&&mod(abs(x1-XE),4)==0&&R<=8 || mod(abs(y1-YE),7)==0&&mod(abs(x1-XE),7)==0&&R>13
                                        %NODE(y1,x1)=1; NGR(y1,x1)=img(y1,x1);
                                    %end
                                else 
                                    sigs = abs(img(y1,x1)-avg);
                                    if sigs < LOCAL_THRESHOLD, LOCAL_THRESHOLD = sigs; end
                                end
                            end
                        end  
                    end
                end
                %
                if id1 > id2
                    EN2=EN2+1;
                    SYe(EN2)=y; SXe(EN2)=x;
                end
            %
            k=k+1;
        end
        %
        if (EN2==0) || (LOCAL_THRESHOLD_MIN>=GLOBAL_THRESHOLD)
            break;
        else
            k=1; EN=EN2;
            if LOCAL_THRESHOLD_MIN>=LOCAL_THRESHOLD
                LOCAL_THRESHOLD_MIN=LOCAL_THRESHOLD_MIN+0.3;
            else
                LOCAL_THRESHOLD_MIN=LOCAL_THRESHOLD+0.1;
            end
            loop=loop+1;
        end
    end
end
%
%
end
%sig0 luon tang len 1 khi bat dau new wave, hoac tang nhieu lan voi 1 den
%khi bat dau sang new wave!
%Thuat toan chuan: Luon co diem dung phan doan khi running neu khong chon all extrema -> khac
%biet so voi thuat toan SRG khac! img(y1,x1)==EXT-sig0*sgn
% -> chi phan doan vung loi (from max) hoac lom (from min) -> full segments if all
% extrema are selected!
%%
% function [Ms, SP, SA, LOCAL_THRESHOLD_MIN, EN2, NODE] = SRG_Nguyen2020_CLERG_full(Ms, SYe, SXe, SG, EN, img)
% % Segmentation algorithm CLERG+LE (proposed by Tuan Nguyen, Tsviatkou, 2020)
% % Name: CLERG+LE full (Model CLERG+LE)
% % Input:
% % - Ms: segmented iamge;
% % - SG: stack of one pixel extreman and their value 
% % - SYe, SXe: stacks of coordinates of the first extreme pixels for Seeded Region Growing
% % - SE: intensity of local extrema
% % - EN: number of extrema
% % - img: input grayscale image
% % Output:
% % - Ms: segmented image of one pixel extrema.
% % - SP: number of loops for speed stack.
% % - LOCAL_THRESHOLD_MIN: threshold
% % - SPEED: speed stack.
% %
% % Convert RGB to grayscale image
% if ndims(img)==3
%     img = rgb2gray(img);
% end
% [Y, X] = size(img);
% img = double(img);
% %
% NODE = sign(abs(Ms)); th=4; SXE=SXe(1:EN); SYE=SYe(1:EN);
% %
% SA=SG;
% %
% B=8; 
% GLOBAL_THRESHOLD = 2^B-1; 
% LOCAL_THRESHOLD_MIN=0; 
% loop=1; 
% k=1;
% EN2=1; SP=0;
% while EN2>0
%     EN2=0; LOCAL_THRESHOLD=2^B-1;
%     while k<=EN
%             y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label); EXT=SG(index,2); sgn=sign(label);
%             id1=0; id2=0;
%             for j=-1:1
%                 for i=-1:1
%                     if (j==0)&&(i==0)
%                     else
%                         y1=y+j; x1=x+i;
%                         if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
%                             id1=id1+1; %avg=SA(index,2);
%                             if img(y1,x1)==EXT-LOCAL_THRESHOLD_MIN*sgn
%                                 Ms(y1,x1)=label; 
%                                 EN=EN+1; SYe(EN)=y1; SXe(EN)=x1;
%                                 idx=SA(index,1);
%                                 %SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
%                                 SA(index,1) = idx+1;
%                                 XE=SXE(index); YE=SYE(index);
% %                                 R=sqrt((y1-YE)^2+(x1-XE)^2);
% %                                 if mod(R,th)==0 && (x1==XE || y1==YE || abs(x1-XE)==abs(y1-YE))
%                                 if mod(abs(y1-YE),th)==0&&mod(abs(x1-XE),th)==0&&(x1==XE || y1==YE || abs(x1-XE)==abs(y1-YE))
%                                     NODE(y1,x1)=1;
%                                 end
%                                 id2=id2+1;
%                             else
%                                 if abs(EXT-img(y1,x1))>LOCAL_THRESHOLD_MIN && abs(EXT-img(y1,x1))<LOCAL_THRESHOLD
%                                     LOCAL_THRESHOLD=abs(EXT-img(y1,x1));
%                                 end   
%                             end
%                         end
%                     end 
%                 end
%             end
%             %
%         if id1>id2
%             EN2=EN2+1;
%             SYe(EN2)=y; SXe(EN2)=x;
%         end
%         %
%         k=k+1;
%     end
%     %
%     if (EN2==0) || (LOCAL_THRESHOLD_MIN>=GLOBAL_THRESHOLD) || (LOCAL_THRESHOLD==2^B-1)
%         break;
%     else
%         k=1; EN=EN2;
%         LOCAL_THRESHOLD_MIN=LOCAL_THRESHOLD;
%         loop=loop+1;
%     end
% end
% %
% k=1;
% while k<=EN2
%     y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label);
%     for j=-1:1
%         for i=-1:1
%             if (j==0)&&(i==0)
%             else
%                 y1=y+j; x1=x+i; 
%                 if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
%                     Ms(y1,x1)=label;
%                     idx=SA(index,1); %avg=SA(index,2);
%                     %SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
%                     SA(index,1) = idx+1;
%                     NODE(y1,x1)=1;%%%%%%%%%%%%%%%%%%%%%
%                     EN2=EN2+1;
%                     SYe(EN2)=y1; SXe(EN2)=x1;
%                 end
%             end
%         end
%     end
%     k=k+1;
% end
% %
% end
% %sig0 luon tang len 1 khi bat dau new wave, hoac tang nhieu lan voi 1 den
% %khi bat dau sang new wave!
% %Thuat toan chuan: Luon co diem dung phan doan khi running neu khong chon all extrema -> khac
% %biet so voi thuat toan SRG khac! img(y1,x1)==EXT-sig0*sgn
% % -> chi phan doan vung loi (from max) hoac lom (from min) -> full segments if all
% % extrema are selected!
%%
% function [Ms, SP, SA, LOCAL_THRESHOLD_MIN, EN2, NODE] = SRG_Nguyen2020_CLERG_full(Ms, SYe, SXe, SG, EN, img)
% % Segmentation algorithm CLERG+LE (proposed by Tuan Nguyen, Tsviatkou, 2020)
% % Name: CLERG+LE full (Model CLERG+LE)
% % Input:
% % - Ms: segmented iamge;
% % - SG: stack of one pixel extreman and their value 
% % - SYe, SXe: stacks of coordinates of the first extreme pixels for Seeded Region Growing
% % - SE: intensity of local extrema
% % - EN: number of extrema
% % - img: input grayscale image
% % Output:
% % - Ms: segmented image of one pixel extrema.
% % - SP: number of loops for speed stack.
% % - LOCAL_THRESHOLD_MIN: threshold
% % - SPEED: speed stack.
% %
% % Convert RGB to grayscale image
% if ndims(img)==3
%     img = rgb2gray(img);
% end
% [Y, X] = size(img);
% img = double(img);
% %
% NODE = sign(abs(Ms)); th=20;
% %
% SA=SG;
% %
% B=8; 
% GLOBAL_THRESHOLD = 2^B-1; 
% LOCAL_THRESHOLD_MIN=0; 
% loop=1; 
% k=1;
% EN2=1; SP=0;
% while EN2>0
%     EN2=0; LOCAL_THRESHOLD=2^B-1;
%     while k<=EN
%             y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label); EXT=SG(index,2); sgn=sign(label);
%             id1=0; id2=0;
%             for j=-1:1
%                 for i=-1:1
%                     if (j==0)&&(i==0)
%                     else
%                         y1=y+j; x1=x+i;
%                         if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
%                             id1=id1+1; %avg=SA(index,2);
%                             if img(y1,x1)==EXT-LOCAL_THRESHOLD_MIN*sgn
%                                 Ms(y1,x1)=label; 
%                                 EN=EN+1; SYe(EN)=y1; SXe(EN)=x1;
%                                 idx=SA(index,1);
%                                 %SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
%                                 SA(index,1) = idx+1;
%                                 if mod((idx+1),th)==0
%                                     NODE(y1,x1)=1;
%                                 end
%                                 id2=id2+1;
%                             else
%                                 if abs(EXT-img(y1,x1))>LOCAL_THRESHOLD_MIN && abs(EXT-img(y1,x1))<LOCAL_THRESHOLD
%                                     LOCAL_THRESHOLD=abs(EXT-img(y1,x1));
%                                 end   
%                             end
%                         end
%                     end 
%                 end
%             end
%             %
%         if id1>id2
%             EN2=EN2+1;
%             SYe(EN2)=y; SXe(EN2)=x;
%         end
%         %
%         k=k+1;
%     end
%     %
%     if (EN2==0) || (LOCAL_THRESHOLD_MIN>=GLOBAL_THRESHOLD) || (LOCAL_THRESHOLD==2^B-1)
%         break;
%     else
%         k=1; EN=EN2;
%         LOCAL_THRESHOLD_MIN=LOCAL_THRESHOLD;
%         loop=loop+1;
%     end
% end
% %
% k=1;
% while k<=EN2
%     y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label);
%     for j=-1:1
%         for i=-1:1
%             if (j==0)&&(i==0)
%             else
%                 y1=y+j; x1=x+i; 
%                 if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
%                     Ms(y1,x1)=label;
%                     idx=SA(index,1); %avg=SA(index,2);
%                     %SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
%                     SA(index,1) = idx+1;
%                     if mod((idx+1),th)==0
%                         NODE(y1,x1)=1;
%                     end
%                     EN2=EN2+1;
%                     SYe(EN2)=y1; SXe(EN2)=x1;
%                 end
%             end
%         end
%     end
%     k=k+1;
% end
% %
% end
% %sig0 luon tang len 1 khi bat dau new wave, hoac tang nhieu lan voi 1 den
% %khi bat dau sang new wave!
% %Thuat toan chuan: Luon co diem dung phan doan khi running neu khong chon all extrema -> khac
% %biet so voi thuat toan SRG khac! img(y1,x1)==EXT-sig0*sgn
% % -> chi phan doan vung loi (from max) hoac lom (from min) -> full segments if all
% % extrema are selected!
%%
% function [Ms, SP, SA, LOCAL_THRESHOLD_MIN, SPEED] = SRG_Nguyen2020_CLERG_full(Ms, SYe, SXe, SG, EN, img)
% % Segmentation algorithm CLERG+LE (proposed by Tuan Nguyen, Tsviatkou, 2020)
% % Name: CLERG+LE full (Model CLERG+LE)
% % Input:
% % - Ms: segmented iamge;
% % - SG: stack of one pixel extreman and their value 
% % - SYe, SXe: stacks of coordinates of the first extreme pixels for Seeded Region Growing
% % - SE: intensity of local extrema
% % - EN: number of extrema
% % - img: input grayscale image
% % Output:
% % - Ms: segmented image of one pixel extrema.
% % - SP: number of loops for speed stack.
% % - LOCAL_THRESHOLD_MIN: threshold
% % - SPEED: speed stack.
% %
% % Convert RGB to grayscale image
% if ndims(img)==3
%     img = rgb2gray(img);
% end
% [Y, X] = size(img);
% img = double(img);
% %
% SA=SG;
% %
% B=8; 
% GLOBAL_THRESHOLD = 2^B-1; 
% LOCAL_THRESHOLD_MIN=0; 
% SPEED=zeros(1,2^B-1);
% loop=1; 
% k=1;
% EN2=1; SP=0;
% while EN2>0
%     EN2=0; LOCAL_THRESHOLD=2^B-1; DeltaE=0; 
%     %DelMAX=0; DelMIN=0; 
%     DeltaN=0; DelNMAX=0; DelNMIN=0;
%     while k<=EN
%             y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label); EXT=SG(index,2); sgn=sign(label);
%             id1=0; id2=0;
%             for j=-1:1
%                 for i=-1:1
%                     if (j==0)&&(i==0)
%                     else
%                         y1=y+j; x1=x+i;
%                         if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
%                             avg=SA(index,2); id1=id1+1;
%                             if img(y1,x1)==EXT-LOCAL_THRESHOLD_MIN*sgn
%                                 Ms(y1,x1)=label; 
%                                 EN=EN+1; SYe(EN)=y1; SXe(EN)=x1;
%                                 idx=SA(index,1);
%                                 SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
%                                 SA(index,1) = idx+1;
%                                 DeltaE = DeltaE+img(y1,x1); 
% %                                 DelMAX=DelMAX+img(y1,x1)*sign(sign(label)+1);
% %                                 DelMIN=DelMIN+img(y1,x1)*sign(1-sign(label));
%                                 DeltaN = DeltaN+1;
%                                 DelNMAX=DelNMAX+sign(sign(label)+1); 
%                                 DelNMIN=DelNMIN+sign(1-sign(label));
%                                 id2=id2+1;
%                             else
%                                 if abs(EXT-img(y1,x1))>LOCAL_THRESHOLD_MIN && abs(EXT-img(y1,x1))<LOCAL_THRESHOLD
%                                     LOCAL_THRESHOLD=abs(EXT-img(y1,x1));
%                                 end   
%                             end
%                         end
%                     end 
%                 end
%             end
%             %
%         if id1>id2
%             EN2=EN2+1;
%             SYe(EN2)=y; SXe(EN2)=x;
%         end
%         %
%         k=k+1;
%     end
%     %
%     if (EN2==0) || (LOCAL_THRESHOLD_MIN>=GLOBAL_THRESHOLD) || (LOCAL_THRESHOLD==2^B-1)
%         break;
%     else
%         k=1; EN=EN2;
%         if DeltaN~=0
%             SP=SP+1;
%             SPEED(SP)=DeltaE/DeltaN;
%         end
%         LOCAL_THRESHOLD_MIN=LOCAL_THRESHOLD;
%         loop=loop+1;
%     end
% end
% %
% SPEED(SP+1:end)=[];
% %
% k=1;
% while k<=EN2
%     y = SYe(k); x = SXe(k); label=Ms(y,x); index=abs(label);
%     for j=-1:1
%         for i=-1:1
%             if (j==0)&&(i==0)
%             else
%                 y1=y+j; x1=x+i; 
%                 if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&Ms(y1,x1)==0
%                     Ms(y1,x1)=label;
%                     idx=SA(index,1); avg=SA(index,2);
%                     SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
%                     SA(index,1) = idx+1;
%                     EN2=EN2+1;
%                     SYe(EN2)=y1; SXe(EN2)=x1;
%                 end
%             end
%         end
%     end
%     k=k+1;
% end
% %
% end