%% Руководство пользователя
% Бессеточная аппроксимация, пространственная фильтрация, текстурная
% сегментация изображений с широким динаническим диапазоном
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% 4.1 БЕССЕТОЧНАЯ АППРОКСИМАЦИЯ ОБЛАСТЕЙ ЛОКАЛЬНЫХ ЭКСТРЕМУМОВ
clear;
A=imread('lena_64x64_8bit.bmp');
%A=imread('Thermal1_64x64_16bit.pgm');
%A=imread('Thermal2_64x64_16bit.pgm');
%A=imread('Thermal3_64x64_16bit.pgm');
%A=imread('Aerial1_64x64_16bit.tif');
%A=imread('Aerial2_64x64_16bit.tif');
%A=imread('Aerial3_64x64_16bit.tif');
%imhist(BD);
%
BD=16; % Битовая глубина
maxx=max(max(A));
if maxx<=255
    BD=8;
end
%
figure, imshow(A, []), impixelinfo;
%% Каскад НЧ-фильтров
img=A;
tic
dim1=3;
im1 = filter_LF_average(img, [dim1 dim1],'symmetric');
dim2=5;
im2 = filter_LF_average(im1, [dim2 dim2],'symmetric');
dim3=7;
im3 = filter_LF_average(im2, [dim3 dim3],'symmetric');
dim4=9;
im4 = filter_LF_average(im3, [dim4 dim4],'symmetric');
toc
figure, imshow(im4, []),title('НЧ информация'),  impixelinfo;
%% Поиск примитивов CSI(4) для первой бессеточной аппроксимации
tic
[SME4, SMM4, EN4, ENs4, SYe4, SXe4, SE4, SG4, OPP4] = Extrema_BSA(im4);
[SMI4, loop4, SA4, sig04, AA4, NODE4, NGR4] = SRG_Nguyen2020_CLERG_appx(SMM4, SYe4, SXe4, SG4, EN4, im4, 10, 20, 15, 29); %(12-24,17-33)
% 
%[SMM4, EN4, SYe4, SXe4, SE4, SG4, OPP4] = Extrema_Single_Symm_Scan(im4);
%[SMI4, loop4, SA4, sig04, AA4, NODE4, NGR4] = SRG_Adams1994_full_quasi(SMM4, SYe4, SXe4, SG4, EN4, im4, 10, 20, 15, 29); %(10-20,15-29)
toc
[r4,c4]=find(NODE4); N4=length(r4) %Число примитивов CSI(4)
figure, imshow(NGR4, []), title('Изображение примитивов'), impixelinfo;
%% Первая бессеточная аппроксимация из CSI(4)
tic
[Y,X]=size(NGR4); REC4=zeros(Y,X);%%%%
M=1; N=1; Dy=ceil(Y/M); Dx=ceil(X/N);
Ny=floor((Y-1)/Dy); Nx=floor((X-1)/Dx);
r=0.4;
for j=0:Ny
    for i=0:Nx
        Ly0=j*Dy+1;             Ly=(j+1)*Dy; 
        Lx0=i*Dx+1;             Lx=(i+1)*Dx;
        ry0=sign(j)*ceil(Dy*r); ry=sign(Ny-j)*ceil(Dy*r);
        rx0=sign(i)*ceil(Dx*r); rx=sign(Nx-i)*ceil(Dx*r);
        if j==Ny, Ly=Y; end
        if i==Nx, Lx=X; end
        %
        [y,x]=find(NODE4(Ly0-ry0:Ly+ry,Lx0-rx0:Lx+rx)); %[rows,cols]%%%%
        y=y+Ly0-ry0-1; x=x+Lx0-rx0-1;
        z=zeros(length(y),1);
        for k=1:length(y)
            z(k)=NGR4(y(k),x(k));%%%%
        end
        %
        [XI,YI] = meshgrid(Lx0:Lx,Ly0:Ly); %(cols, rows)
        ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'invquadratic', 'RBFConstant', 10)); %12, 10, 8, 6, 4
        ZI = reshape(ZI, size(XI));
        REC4(Ly0:Ly,Lx0:Lx)=ZI;%%%%
        %
    end
end
toc
[mse, rmse, psnr] = ipsnr(im4, REC4, BD)%%%%
figure, imshow(REC4, []), title('Восстановленное изображение'), impixelinfo;%%%%
%% Поиск примитивов CSI(3) для второй бессеточной аппроксимации
E3 = im3-REC4; figure, imshow(E3, []),  title('Ошибка аппроксимации'), impixelinfo;
[SME3, SMM3, EN3, ENs3, SYe3, SXe3, SE3, SG3, OPP3] = Extrema_BSA(E3);
[SMI3, loop3, SA3, sig03, AA3, NODE3, NGR3] = SRG_Nguyen2020_CLERG_appx(SMM3, SYe3, SXe3, SG3, EN3, E3, 8, 16, 13, 25); %(8-16,13-25)
%
%[SMM3, EN3, SYe3, SXe3, SE3, SG3, OPP3] = Extrema_Single_Symm_Scan(E3);
%[SMI3, loop3, SA3, sig03, AA3, NODE3, NGR3] = SRG_Adams1994_full_quasi(SMM3, SYe3, SXe3, SG3, EN3, E3, 8, 16, 13, 25); %(8-16,13-25)
[r3,c3]=find(NODE3); N3=length(r3) %Число примитивов CSI(3)
figure, imshow(NODE3, []), title('Изображение примитивов'), impixelinfo;
%% Вторая бессеточная аппроксимация из CSI(3)
tic
[Y,X]=size(NGR3); REC3=zeros(Y,X);%%%%
M=1; N=1; Dy=ceil(Y/M); Dx=ceil(X/N);
Ny=floor((Y-1)/Dy); Nx=floor((X-1)/Dx);
r=0.4;
for j=0:Ny
    for i=0:Nx
        Ly0=j*Dy+1;             Ly=(j+1)*Dy; 
        Lx0=i*Dx+1;             Lx=(i+1)*Dx;
        ry0=sign(j)*ceil(Dy*r); ry=sign(Ny-j)*ceil(Dy*r);
        rx0=sign(i)*ceil(Dx*r); rx=sign(Nx-i)*ceil(Dx*r);
        if j==Ny, Ly=Y; end
        if i==Nx, Lx=X; end
        %
        [y,x]=find(NODE3(Ly0-ry0:Ly+ry,Lx0-rx0:Lx+rx)); %[rows,cols]%%%%
        y=y+Ly0-ry0-1; x=x+Lx0-rx0-1;
        z=zeros(length(y),1);
        for k=1:length(y)
            z(k)=NGR3(y(k),x(k));%%%%
        end
        %
        [XI,YI] = meshgrid(Lx0:Lx,Ly0:Ly); %(cols, rows)
        ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'invquadratic', 'RBFConstant', 6));
        ZI = reshape(ZI, size(XI));
        REC3(Ly0:Ly,Lx0:Lx)=ZI;%%%%
        %
    end
end
toc
[mse, rmse, psnr] = ipsnr(E3, REC3, BD)%%%%
figure, imshow(REC3, []), title('Восстановленное изображение'), impixelinfo;%%%%
%% Поиск примитивов CSI(2) для третьей бессеточной аппроксимации
E2=im2-REC4-REC3; figure, imshow(E2, []), title('Ошибка аппроксимации'), impixelinfo;
[SME2, SMM2, EN2, ENs2, SYe2, SXe2, SE2, SG2, OPP2] = Extrema_BSA(E2);
[SMI2, loop2, SA2, sig02, AA2, NODE2, NGR2] = SRG_Nguyen2020_CLERG_appx(SMM2, SYe2, SXe2, SG2, EN2, E2, 6, 12, 9, 17);%(6-12,9-17)
%
%[SMM2, EN2, SYe2, SXe2, SE2, SG2, OPP2] = Extrema_Single_Symm_Scan(E2);
%[SMI2, loop2, SA2, sig02, AA2, NODE2, NGR2] = SRG_Adams1994_full_quasi(SMM2, SYe2, SXe2, SG2, EN2, E2, 6, 12, 9, 17); %(6-12,9-17)
[r2,c2]=find(NODE2); N2=length(r2) %Число примитивов CSI(2)
figure, imshow(NODE2, []), title('Изображение примитивов'), impixelinfo;
%% Третья бессеточная аппроксимация из CSI(2)
tic
[Y,X]=size(NGR2); REC2=zeros(Y,X);%%%%
M=1; N=1; Dy=ceil(Y/M); Dx=ceil(X/N);
Ny=floor((Y-1)/Dy); Nx=floor((X-1)/Dx);
r=0.4;
for j=0:Ny
    for i=0:Nx
        Ly0=j*Dy+1;             Ly=(j+1)*Dy; 
        Lx0=i*Dx+1;             Lx=(i+1)*Dx;
        ry0=sign(j)*ceil(Dy*r); ry=sign(Ny-j)*ceil(Dy*r);
        rx0=sign(i)*ceil(Dx*r); rx=sign(Nx-i)*ceil(Dx*r);
        if j==Ny, Ly=Y; end
        if i==Nx, Lx=X; end
        %
        [y,x]=find(NODE2(Ly0-ry0:Ly+ry,Lx0-rx0:Lx+rx)); %[rows,cols]%%%%
        y=y+Ly0-ry0-1; x=x+Lx0-rx0-1;
        z=zeros(length(y),1);
        for k=1:length(y)
            z(k)=NGR2(y(k),x(k));%%%%
        end
        %
        [XI,YI] = meshgrid(Lx0:Lx,Ly0:Ly); %(cols, rows)
        ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'invquadratic', 'RBFConstant', 4));
        ZI = reshape(ZI, size(XI));
        REC2(Ly0:Ly,Lx0:Lx)=ZI;%%%%
        %
    end
end
toc
[mse, rmse, psnr] = ipsnr(E2, REC2, BD)%%%%
figure, imshow(REC2, []), title('Восстановленное изображение'), impixelinfo;%%%%
%% Поиск примитивов CSI(1) для четворой бессеточной аппроксимации
E1=im1-REC4-REC3-REC2; figure, imshow(E1, []), title('Ошибка аппроксимации'), impixelinfo;
[SME1, SMM1, EN1, ENs1, SYe1, SXe1, SE1, SG1, OPP1] = Extrema_BSA(E1);
[SMI1, loop1, SA1, sig01, AA1, NODE1, NGR1] = SRG_Nguyen2020_CLERG_appx(SMM1, SYe1, SXe1, SG1, EN1, E1, 4, 8, 7, 13); %(4-8,7-13)
%
%[SMM1, EN1, SYe1, SXe1, SE1, SG1, OPP1] = Extrema_Single_Symm_Scan(E1);
%[SMI1, loop1, SA1, sig01, AA1, NODE1, NGR1] = SRG_Adams1994_full_quasi(SMM1, SYe1, SXe1, SG1, EN1, E1, 4, 8, 7, 13); %(4-8,7-13)
[r1,c1]=find(NODE1); N1=length(r1) %Число примитивов CSI(1)
figure, imshow(NODE1, []), title('Изображение примитивов'), impixelinfo;
%% Четвортая бессеточная аппроксимация из CSI(1)
tic
[Y,X]=size(NGR1); REC1=zeros(Y,X);%%%%
M=1; N=1; Dy=ceil(Y/M); Dx=ceil(X/N);
Ny=floor((Y-1)/Dy); Nx=floor((X-1)/Dx);
r=0.4;
for j=0:Ny
    for i=0:Nx
        Ly0=j*Dy+1;             Ly=(j+1)*Dy; 
        Lx0=i*Dx+1;             Lx=(i+1)*Dx;
        ry0=sign(j)*ceil(Dy*r); ry=sign(Ny-j)*ceil(Dy*r);
        rx0=sign(i)*ceil(Dx*r); rx=sign(Nx-i)*ceil(Dx*r);
        if j==Ny, Ly=Y; end
        if i==Nx, Lx=X; end
        %
        [y,x]=find(NODE1(Ly0-ry0:Ly+ry,Lx0-rx0:Lx+rx)); %[rows,cols]%%%%
        y=y+Ly0-ry0-1; x=x+Lx0-rx0-1;
        z=zeros(length(y),1);
        for k=1:length(y)
            z(k)=NGR1(y(k),x(k));%%%%
        end
        %
        [XI,YI] = meshgrid(Lx0:Lx,Ly0:Ly); %(cols, rows)
        ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'invquadratic', 'RBFConstant', 3));
        ZI = reshape(ZI, size(XI));
        REC1(Ly0:Ly,Lx0:Lx)=ZI;%%%%
        %
    end
end
toc
[mse, rmse, psnr] = ipsnr(E1, REC1, BD)%%%%
figure, imshow(REC1, []), title('Восстановленное изображение'), impixelinfo;%%%%
%% Ошибка аппроксимации
E0=double(img)-REC4-REC3-REC2-REC1; figure, imshow(E0, []), title('Ошибка аппроксимации E0'), impixelinfo;
%% Восстановленное изображение
figure, imshow(REC4+REC3+REC2+REC1, []), title('Восстановленное изображение A1'), impixelinfo;
%%
figure, imshow(REC4+REC3+REC2, []), title('Восстановленное изображение A2'), impixelinfo;
%%
figure, imshow(REC4+REC3, []), title('Восстановленное изображение A3'), impixelinfo;
%%
figure, imshow(REC4, []), title('Восстановленное изображение A4'), impixelinfo;
%% Исходное изображение
A00 = REC4+REC3+REC2+REC1+E0; figure, imshow(A00, []), title('Исходное изображение A2'), impixelinfo;
%% Энергия аппроксимации
sum(sum(abs(NGR4)))+ sum(sum(abs(NGR3)))+sum(sum(abs(NGR2)))+sum(sum(abs(NGR1)))+sum(sum(abs(E0)))
%% Энергия
sum(sum(abs(NGR4))) % 4-ый уровень
%% Энергия
sum(sum(abs(NGR3))) % 3-й уровень
%% Энергия
sum(sum(abs(NGR2))) % 2-ой уровень
%% Энергия
sum(sum(abs(NGR1))) % 1-ый уровень
%% Энергия
sum(sum(abs(E0)))   % Ошибка аппроксимации
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%% 4.2 Пространственная фильтрация импульсных шумов и малоразмерных объектов на изображениях
% 4.2.1. Удаление импульсных шумов на основе анализа яркости областей ЛЭ
%
% 4.2.2. Удаление малоразмерных объектов на основе анализа размеров областей ЛЭ
%--------------------------------------------------------------------------
%% 4.2.1. Удаление импульсных шумов на основе анализа яркости областей ЛЭ
%clear;
%I=imread('lena_256x256_8bit.bmp');
I=imread('Thermal1_320x256_16bit.pgm');
%I=imread('Thermal2_320x256_16bit.pgm');
%I=imread('Thermal3_320x256_16bit.pgm');
%I=imread('Aerial1_320x256_16bit.tif');
%I=imread('Aerial2_320x256_16bit.tif');
%I=imread('Aerial3_320x256_16bit.tif');
BD=16; % Битовая глубина
maxx=max(max(I));
if maxx<=255
    BD=8;
end
%% Зашумление изображений при различной плотности импульсного шума
tic
J = imnoise(I,'salt & pepper',0.2); %0.02 - noise density = d*numel(I)
if BD==16
    minn = min(min(I)); maxx=max(max(I));
    J(J==0)=minn; J(J==65535)=maxx; %Для 16-битных изображений
else
    minn = 0; maxx=255; %Для 8-битных изображений
end
toc
%
%figure, imshow(I, []), title('Исходное изображение'), impixelinfo;
%figure, imshow(J, []), title('Зашумленное изображение'), impixelinfo;
figure, imshowpair(I, J, 'montage'), title('Исходное - Зашумленное'), impixelinfo;
%% Удаление импульсного шума Медианным фильтром
tic
K = medfilt2(J, [3 3]);
K(K==0)=minn;
toc
%
[mse, rmse, psnr] = ipsnr(I, K, BD) %BD= 16bit, 8bit
%figure, imshow(K, []), title('Восстановленное изображение'), impixelinfo;
figure, imshowpair(J, K, 'montage'), title('Зашумленное - Восстановленное'), impixelinfo;
%% Удаление импульсного шума на основе блочного поиска экстремумов
tic
[SME, EN, SYe, SXe, SE, SG, OPP] = Extrema_Single_Symm_Scan(J);
[imF, imS, SME2, NSS, SYs, SXs] = filter_salt_pepper(J, SME, SYe, SXe, SE, EN, minn, maxx); %16bit: (0,65535); 8bit: (0,255)
toc
%
[mse, rmse, psnr] = ipsnr(I, imF, BD) %BD= 16bit, 8bit
%figure, imshow(imF, []), title('Восстановленное изображение'), impixelinfo;
figure, imshowpair(J, imF, 'montage'), title('Зашумленное - Восстановленное'), impixelinfo;
%% Удаление импульсного шума на основе блочно-сегментного поиска экстремумов
tic
[SME, SMM, EN, ENs, SYe, SXe, SE, SG, OPP] = Extrema_BSA(J);
[imF, imS, SME2, NSS, SYs, SXs] = filter_salt_pepper(J, SME, SYe, SXe, SE, EN, minn, maxx); %16bit: (0,65535); 8bit: (0,255)
toc
%
[mse, rmse, psnr] = ipsnr(I, imF, BD) %BD= 16bit, 8bit
%figure, imshow(imF, []), title('Восстановленное изображение'), impixelinfo;
figure, imshowpair(J, imF, 'montage'), title('Зашумленное - Восстановленное'), impixelinfo;
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%% 4.2.2. Удаление малоразмерных объектов на основе анализа размеров областей ЛЭ
%clear;
img=imread('Thermal1_320x256_16bit.pgm');
%img=imread('Thermal2_320x256_16bit.pgm');
%img=imread('Thermal3_320x256_16bit.pgm');
%img=imread('Aerial1_320x256_16bit.tif');
%img=imread('Aerial2_320x256_16bit.tif');
%img=imread('Aerial3_320x256_16bit.tif');
%img=imread('Aerial4_256x256_8bit.tif');
BD=16; % Битовая глубина
maxx=max(max(img));
if maxx<=255
    BD=8;
end
figure, imshow(img, []), title('Исходное'), impixelinfo;
% Вычисление среднего градиента исходного изображения
I=img;
[Gx,Gy] = imgradientxy(I);
[Gmag,Gdir] = imgradient(Gx,Gy);
AG = sum(sum(Gmag))./(sqrt(2)*(size(I,1)-1)*(size(I,2)-1))
%% Удаление малоразмерных объектов на основе блочного поиска экстремумов
tic
[SMM, EN, SYe, SXe, SE, SG, OPP] = Extrema_Single_Symm_Scan(img);
[SMI, loop, SA, sig0, AA] = SRG_Nguyen2020_CLERG_object(SMM, SYe, SXe, SG, EN, img);
SIZE=10; %Размер обектов
[imF, imS, SMI2, NSS, SYs, SXs] = filter_object(img, SMI, SYe, SXe, SA, SE, EN, SIZE, maxx, minn);
toc
[mse, rmse, psnr] = ipsnr(img, imF, BD) %16bit: 65535; 8bit: 255
%figure, imshow(imS, []), title('Малоразмерные объекты'), impixelinfo;
figure, imshow(imF, []), title('Восстановленное'), impixelinfo;
%
% Вычисление среднего градиента
I=imF;
[Gx,Gy] = imgradientxy(I);
[Gmag,Gdir] = imgradient(Gx,Gy);
AG = sum(sum(Gmag))./(sqrt(2)*(size(I,1)-1)*(size(I,2)-1))
%% Удаление малоразмерных объектов на основе блочно-сегментного поиска экстремумов
tic
[SME, SMM, EN, ENs, SYe, SXe, SE, SG, OPP] = Extrema_BSA(img);
[SMI, loop, SA, sig0, AA] = SRG_Nguyen2020_CLERG_object(SMM, SYe, SXe, SG, EN, img);
SIZE=10; %Размер обектов
[imF, imS, SMI2, NSS, SYs, SXs] = filter_object(img, SMI, SYe, SXe, SA, SE, EN, SIZE, maxx, minn);
toc
[mse, rmse, psnr] = ipsnr(img, imF, BD) %16bit: 65535; 8bit: 255
%figure, imshow(imS, []), title('Малоразмерные объекты'), impixelinfo;
figure, imshow(imF, []), title('Восстановленное'), impixelinfo;
%
% Вычисление среднего градиента
I=imF;
[Gx,Gy] = imgradientxy(I);
[Gmag,Gdir] = imgradient(Gx,Gy);
AG = sum(sum(Gmag))./(sqrt(2)*(size(I,1)-1)*(size(I,2)-1))
%% Удаление малоразмерных объектов Гауссовским фильтром
I=img;
h = fspecial('gaussian', [11 11], 5);
tic
imG = imfilter(I, h, 'symmetric');
toc
[mse, rmse, psnr] = ipsnr(I, imG, BD) %16bit: 65535; 8bit: 255
figure, imshow(imG, []), title('Отфильтрованное'), impixelinfo;
%
% Вычисление среднего градиента
I=imG;
[Gx,Gy] = imgradientxy(I);
[Gmag,Gdir] = imgradient(Gx,Gy);
AG = sum(sum(Gmag))./(sqrt(2)*(size(I,1)-1)*(size(I,2)-1))
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%% 4.2.3 Текстурная сегментация на основе оценки плотности областей ЛЭ
% 1. Для 8-битных синтетических текстурных изображений
img=imread('01.png');
%img=imread('02.png');
%img=imread('03.png');
%img=imread('04.png');
%img=imread('07.bmp'); %img = rgb2gray(img); immm=img(1:end-1,1:end-1); imwrite(immm, '07.bmp');
%img=imread('06.png');
if ndims(img)==3
    img = rgb2gray(img);
end
tic
[SME, SMM, EN, ENs, SYe, SXe, SE, SG, OPP] = Extrema_BSA(img);
%[SME, EN, SYe, SXe, SE, SG, OPP] = Extrema_Single_Symm_Scan(img);
SIZE = 60; %Размер окна = 41, 45, 55, 65, 85, 105, 115, 125
%[imINS, imS, SMEE] = filter_extreme_energy(img, SME, SIZE);
[imINS, imS, SMEE] = filter_extreme_energy(img, SME, SIZE);
toc
figure, imshow(img, []), title('Исходное'), impixelinfo;
%figure, imshow(imINS, []), title('Плотность энергии ЛЭ'), impixelinfo;
figure, imshow(imS, []), title('Сегментированное'), impixelinfo;
%% 2. Для 16-битных текстурных изображений
immm = imread('07.tif');
figure, imshow(immm, []),  impixelinfo;
% Текстурная сегментация на основе блочного и блочно-сегментного поиска ЛЭ
tic
[SME, SMM, EN, ENs, SYe, SXe, SE, SG, OPP] = Extrema_BSA(immm);
%[SME, EN, SYe, SXe, SE, SG, OPP] = Extrema_Single_Symm_Scan(immm);
SIZE = 60; %Размер окна = 41, 45, 55, 65, 85, 105, 115, 125
[imINS,imS, SMEE] = filter_extreme_energy(immm, SME, SIZE);
toc
%figure, imshow(imINS, []), title('Плотность энергии ЛЭ'), impixelinfo;
figure, imshow(imS, []), title('Сегментированное'), impixelinfo;
%%