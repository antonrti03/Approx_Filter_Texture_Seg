%https://ch.mathworks.com/help/matlab/math/interpolating-scattered-data.html
%https://ch.mathworks.com/help/matlab/ref/scatteredinterpolant.html#d122e1177219
clear;
X = [-1.5 3.2; 1.8 3.3; -3.7 1.5; -1.5 1.3; ...
      0.8 1.2; 3.3 1.5; -4.0 -1.0; -2.3 -0.7; 
      0 -0.5; 2.0 -1.5; 3.7 -0.8; -3.5 -2.9; ...
      -0.9 -3.9; 2.0 -3.5; 3.5 -2.25; ...
      -4 4; 4 4; -4 -4; 4 -4];
V = X(:,1).^2 + X(:,2).^2; %X(L,1)-Hoanh, X(:,2)-Tung
hold on
plot3(X(:,1),X(:,2),zeros(19,1), '*r')
axis([-4, 4, -4, 4, 0, 32]);
stem3(X(:,1),X(:,2),V,'^','fill')
grid
hold off
view(322.5, 30);
%% Interpolate scattered data
x = X(:,1);
y = X(:,2);
v = V;
[xq,yq] = meshgrid(-4:.1:4, -4:.1:4);
vq = griddata(x,y,v,xq,yq,'natural');
figure
mesh(xq,yq,vq);
hold on
plot3(x,y,v,'or');
%% Create a 2-D or 3-D triangulation from a set of points
DT = delaunayTriangulation(X);
figure
triplot(DT, 'k')
%% Create a Delaunay triangulation, lift the vertices, and evaluate the interpolant at the query point Xq.
figure('Color', 'white')
t = delaunay(X(:,1),X(:,2));
%
hold on
trimesh(t,X(:,1),X(:,2), zeros(15,1), ...
    'EdgeColor','r', 'FaceColor','none')
defaultFaceColor  = [0.6875 0.8750 0.8984];
trisurf(t,X(:,1),X(:,2), V, 'FaceColor', ...
    defaultFaceColor, 'FaceAlpha',0.9);
plot3(X(:,1),X(:,2),zeros(15,1), '*r')
axis([-4, 4, -4, 4, 0, 25]);
grid
plot3(-2.6,-2.6,0,'*b','LineWidth', 1.6)
plot3([-2.6 -2.6]',[-2.6 -2.6]',[0 13.52]','-b','LineWidth',1.6)
hold off
%
view(322.5, 30);
text(-2.0, -2.6, 'Xq', 'FontWeight', 'bold', ...
'HorizontalAlignment','center', 'BackgroundColor', 'none');
%% Interolate Scattered Data Over a Uniform Grid
xy = -2.5 + 5*gallery('uniformdata',[200 2],0);
x = xy(:,1);
y = xy(:,2);
v = x.*exp(-x.^2-y.^2);
[xq,yq] = meshgrid(-2:.2:2, -2:.2:2);
vq = griddata(x,y,v,xq,yq);
figure
mesh(xq,yq,vq);
hold on
plot3(x,y,v,'o');

h = gca;
h.XLim = [-2.7 2.7];
h.YLim = [-2.7 2.7];
%%












%% Interpolate scattered data: Matlab griddata interpolation
x = E(1:NE,2); %colums
y = E(1:NE,1); %rows
v = E(1:NE,3); %values
[xq,yq] = meshgrid(1:1:256, 1:1:256);
tic
% 'linear' (df), 'cubic', 'multiquadric', 'thinplate', 'gaussian'
vq = griddata(x,y,v,xq,yq,'linear');
toc
figure
mesh(yq,xq,vq);
hold on
plot3(y,x,v,'.r');
view(50, 60);
%%
x = E(1:NE,2); %colums
y = E(1:NE,1); %rows
v = E(1:NE,3); %values
[xq,yq] = meshgrid(1:1:256, 1:1:256);
%figure, mesh(xq,yq,im2.');
figure, surf(xq,yq,im2.', 'edgecolor', 'none')
hold on
plot3(y,x,v,'.r');
xlabel('rows'); ylabel('columns');
title('Input image and red nodes');
view(75,45);
%%
vq1=vq(5:256,5:256);
%% Create a 2-D or 3-D triangulation from a set of points
X = [E(:,2) E(:,1)];
DT = delaunayTriangulation(X);
figure
triplot(DT, 'k')
%%




%% RBF Interpolation and Approximation
% G:\PhD1721\Matlab_Approximation\RBF
%https://ch.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions
x = rand(50,1)*4-2; 
y = rand(50,1)*4-2;
z = x.*exp(-x.^2-y.^2);

ti = -2:.05:2; 
[XI,YI] = meshgrid(ti,ti);
%%
ZIg = griddata(x,y,z,XI,YI,'cubic');

%RBF interpolation
ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'multiquadric', 'RBFConstant', 1));

ZI = reshape(ZI, size(XI));

%Plot data
subplot(2,2,1); mesh(XI,YI,ZIg), hold, axis([-2 2 -2 2 -0.5 0.5]); 
plot3(x,y,z,'.r'), hold off; title('Interpolation using Matlab function griddata(method=cubic)');

subplot(2,2,3); pcolor(abs(ZIg - XI.*exp(-XI.^2-YI.^2))); colorbar; title('griddata(method=cubic) interpolation error');

subplot(2,2,2); mesh(XI,YI,ZI), hold
plot3(x,y,z,'.r'), hold off; title('RBF interpolation'); axis([-2 2 -2 2 -0.5 0.5]);
subplot(2,2,4); pcolor(abs(ZI - XI.*exp(-XI.^2-YI.^2))); colorbar; title('RBF interpolation error');

%%
x = E(1:NE,2); 
y = E(1:NE,1); 
z = E(1:NE,3);

ti = 1:1:256; 
[XI,YI] = meshgrid(ti,ti);
% 'linear', 'natural', 'cubic', 'v4'
ZIg = griddata(x,y,z,XI,YI,'cubic');
% Weight values
%coeff = rbfcreate([x'; y'], z','RBFFunction', 'multiquadric','RBFConstant', 25);
%RBF interpolation
tic
% 'linear' (df), 'cubic', 'multiquadric', 'thinplate', 'gaussian'
ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'multiquadric', 'RBFConstant', 25));
% Optimal value of the parameter ?  is somewhat close to the average distance between interpolation nodes.
ZI = reshape(ZI, size(XI));
toc
%Plot data
figure, mesh(XI,YI,ZIg), hold, axis([0 256 0 256 0 255]); 
plot3(x,y,z,'.r'), hold off; title('Interpolation using Matlab function griddata(method=cubic)');

figure, pcolor((im2-ZIg)); colorbar; title('griddata(method=cubic) interpolation error');

figure, mesh(XI,YI,ZI), hold
plot3(x,y,z,'.r'), hold off; title('RBF interpolation'); axis([0 256 0 256 0 255]);
figure, pcolor((im2-ZI)); colorbar; title('RBF interpolation error');
%%
x = E(1:NE,2); %colums
y = E(1:NE,1); %rows
z = E(1:NE,3);
ti = 1:1:256; 
[XI,YI] = meshgrid(ti,ti);
% 'linear', 'natural', 'cubic', 'v4'
%ZIg = griddata(x,y,z,XI,YI,'cubic');
%RBF interpolation
tic
% 'linear' (=df), 'cubic', 'multiquadric', 'thinplate', 'gaussian'
% 'invquadratic', 'invmultiquadric'
ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'invquadratic', 'RBFConstant', 25));
%coeff = rbfcreate([x'; y'], z','RBFFunction', 'invquadratic','RBFConstant', 32);
% Optimal value of the parameter ?  is somewhat close to the average distance between interpolation nodes.
ZI = reshape(ZI, size(XI));
toc
%%
tic
NODE1=NODE(:,37:end);
[y,x]=find(NODE1); %[rows,cols]
x=x+36;
z=zeros(length(y),1);
for k=1:length(y)
    z(k)=img(y(k),x(k));
end
toc
% [Y,X]=size(img);
% y=zeros(Y*X,1); x=zeros(Y*X,1); z=zeros(Y*X,1);
% k=0;
% tic
% for j=1:Y
%     for i=1:X
%         if NODE(j,i)==1
%             k=k+1;
%             y(k)=j; x(k)=i; z(k)=img(j,i);
%         end
%     end
% end
% toc
[XI,YI] = meshgrid(1:1:72, 1:1:48); %(cols, rows)
% 'linear', 'natural', 'cubic', 'v4'
%ZIg = griddata(x,y,z,XI,YI,'cubic');
%RBF interpolation
%%
tic
% 'linear' (=df), 'cubic', 'multiquadric', 'thinplate', 'gaussian'
% 'invquadratic', 'invmultiquadric'
ZI = rbfinterp([XI(:)'; YI(:)'], rbfcreate([x'; y'], z','RBFFunction', 'invquadratic', 'RBFConstant', 7));
%coeff = rbfcreate([x'; y'], z','RBFFunction', 'invquadratic','RBFConstant', 32);
% Optimal value of the parameter ?  is somewhat close to the average distance between interpolation nodes.
ZI = reshape(ZI, size(XI));
toc
%%
%figure, mesh(xq,yq,im2.');
figure, surf(XI, YI, ZI.', 'edgecolor', 'none')
hold on
plot3(y,x,z,'.r');
xlabel('rows'); ylabel('columns');
title('Approx. image and red nodes');
view(75,45);
%% RBF plot
x = 0:0.01:3; 
ym = sqrt(x.*x+1);
yg = exp(-x.*x/2);
yc = x.*x.*x;
yth = x.*x.*log(x+1);
KF = 1;
yp = 1./(1+KF*(x.*x));
plot(x, yp, x, yc, x, yth, x,ym, x, yg, x, x); legend('proposed','cubic', 'thinplate', 'multiquadrics', 'gaussian', 'linear');
grid on
%% RBF plot
x = -2:0.01:2; 
ym = sqrt(x.*x+1);
yim = 1./sqrt(x.*x+1);
yg = exp(-x.*x/2);
KF = 1;
yp = 1./(1+KF*(x.*x)); %proposed
plot(x, yim, x, yp, x,ym, x, yg, x, x); legend('invmultiquadrics', 'proposed', 'multiquadrics', 'gaussian', 'linear');
grid on
%%