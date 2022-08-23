function [Ms, loop, SA, LOCAL_THRESHOLD_MIN, EN2, NODE, NGR] = SRG_Adams1994_full_quasi(Ms, SYe, SXe, SG, EN, img, a, b, c, d)
% Segmentation algorithm: SRG_Adams1994_full_quasi
% Made by Tuan Nguyen, Tsvetkov, 2019
% Name: SRG_Adams1994_full_quasi+LE
% Input:
% - Ms: segmented iamge;
% - SG: stack of one pixel extreman and their value 
% - SYe, SXe: stacks of coordinates of the first extreme pixels for Seeded Region Growing
% - SE: intensity of local extrema
% - EN: number of extrema
% - img: input grayscale image
% Output:
% - Ms: segmented image of one pixel extrema.
% - loop: loop number of SRG.
% - LOCAL_THRESHOLD_MIN: threshold
% - SA: average values of segments
% - SPEED: Check stack.
%
% Rolf Adams, Leanne Bischof. Seeded region growing. IEEE Transactions on Pattern Analysis and Machine Intelligence. 1994; 16(6):641–647.
%
% Convert RGB to grayscale image
if ndims(img)==3
    img = rgb2gray(img);
end
[Y, X] = size(img);
%img = double(img);
img = round(double(img));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
NODE = sign(abs(Ms)); SXE=SXe; SYE=SYe;
NODE(1,1)=1; NODE(end:end)=1; NODE(1,end)=1; NODE(end,1)=1;
NGR  = img.*NODE;
%
SA=SG;
%
maxx=max(max(double(img)));%B=8; 
GLOBAL_THRESHOLD = maxx;
LOCAL_THRESHOLD_MIN=0; 
loop=1; 
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
                            %if abs(img(y1,x1)-img(y,x)) <= LOCAL_THRESHOLD_MIN
                            if abs(img(y1,x1)-avg) <= LOCAL_THRESHOLD_MIN    
                                Ms(y1,x1)=label; 
                                EN=EN+1; SYe(EN)=y1; SXe(EN)=x1;
                                idx=SA(index,1); 
                                SA(index,2) = (avg*idx + img(y1,x1))/(idx+1);
                                SA(index,1) = idx+1;
                                id2=id2+1;
                                XE=SXE(index); YE=SYE(index);
                                R = round(sqrt((x1-XE)^2+(y1-YE)^2));
                                %(10-20,15-29)
                                %(8-16,13-25)
                                %(6-12,9-17)
                                %(4-8,7-13)
                                if mod(abs(y1-YE),a)==0&&mod(abs(x1-XE),a)==0&&R<=b || mod(abs(y1-YE),c)==0&&mod(abs(x1-XE),c)==0&&R>d
                                    NODE(y1,x1)=1; NGR(y1,x1)=img(y1,x1);
                                end
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
% Sig0 changes value to THRESHOLD
%Full segmentation: khong phu thuoc vao initial points!