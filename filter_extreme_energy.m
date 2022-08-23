function [imINS,imS, SMEE] = filter_extreme_energy(img, SME, SIZE)
if ndims(img)==3
    img = rgb2gray(img);
end
img = double(img); img(img==0)=1;
SME = abs(sign(SME)); SMEE=SME.*img;
PAD=floor(SIZE/2);
[Y,X]=size(SMEE);
imINS=zeros(Y,X);
for y=1:Y
    for x=1:X
        yp1=y-PAD; xp1=x-PAD; yp2=y+PAD; xp2=x+PAD;
        if yp1<1, yp1=1; end
        if yp2>Y, yp2=Y; end
        if xp1<1, xp1=1; end
        if xp2>X, xp2=X; end
        imINS(y,x)=sum(sum(SMEE(yp1:yp2, xp1:xp2)))/((yp2-yp1+1)*(xp2-xp1+1));
    end
end
%Thresholding
imnor = imINS./(max(max(imINS))+eps);
imS = im2bw(imnor, graythresh(imnor));
%
end