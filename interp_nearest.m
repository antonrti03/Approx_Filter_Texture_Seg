function [IM_nearest] = interp_nearest(img, NODE)
% Created by Tuan Nguyen, 02/04/2021\
%
[row,col]=find(NODE==1); [Y,X] =size(img);
EY = zeros(Y*X,1); EX = zeros(Y*X,1); 
EY(1:length(row))= row; EX(1:length(col))= col;
IM_nearest = NODE.*double(img);
ENR = length(col);
imcheck = NODE;
k=1;
while k<=ENR
    y = EY(k); x = EX(k);
    for j=-1:1
        for i=-1:1
            if (j==0)&&(i==0)
            else
                y1=y+j; x1=x+i;
                if (y1>0)&&(x1>0)&&(y1<=Y)&&(x1<=X)&&imcheck(y1,x1)==0
                    IM_nearest(y1,x1)=IM_nearest(y,x);
                    imcheck(y1,x1)=1;
                    ENR=ENR+1; EY(ENR)=y1; EX(ENR)=x1;
                end
            end
        end
    end
    k=k+1;
end
%
end