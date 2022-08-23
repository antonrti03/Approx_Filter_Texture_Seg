Radial base functions (RBF) can be used for interpolation and and approximation of scattered data i.e. data is not required to be on any regular grid. The same function can handle data interpolation in any dimension. See file rbftest.m for more examples.

1. Create RBF interpolation using
rbf=rbfcreate(x, f); ?x? ? coordinates of the nodes and ?f? - values of the function at the nodes

2. Calculate interpolated values ?fi? at nodes ?xi? using
fi = rbfinterp(xi, rbf); rbf ? is structure returned by rbf=rbfcreate(x, f)

%1D example
x = 0:1.25:10; f = sin(x);
xi = 0:.1:10;

%Matlab interpolation
fi = interp1(x,f,xi);

% RBF interpolation
rbf=rbfcreate(x, f);
fi = rbfinterp(xi, rbf);

%2D example
x = rand(50,1)*4-2; y = rand(50,1)*4-2; z = x.*exp(-x.^2-y.^2);

ti = -2:.05:2;
[XI,YI] = meshgrid(ti,ti);

%Matlab interpolation
ZI = griddata(x,y,z,XI,YI,'cubic');

%RBF interpolation
rbf=rbfcreate([x'; y'], z');
ZI = rbfinterp([XI(:)'; YI(:)'], op);
ZI = reshape(ZI, size(XI));

Optional parameters:

1. Radial Base Function:
rbfcreate(x, f ,'RBFFunction', 'multiquadric');
available RBF functions are: multiquadric, gaussian, linear, cubic, thinplate
2. Smoothing level: (must be a positive scalar)
rbfcreate(x, f ,'RBFSmooth', 0.1);
3. Multiquadric and gaussian functions have definable constants
rbfcreate(x, f ,?RBFConstant', 0.1);

RBF interpolation usually produces much better results that standard Matlab functions but computation complexity of RBF interpolation is n^3 thus it is not recommended to use it for more then 2000 nodes.

Near-optimal data-independent point locations for radial basis function interpolation
================
%% RBF Interpolation and Approximation
% G:\PhD1721\Matlab_Approximation\RBF
%https://ch.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions
x = rand(50,1)*4-2; 
y = rand(50,1)*4-2;
z = x.*exp(-x.^2-y.^2);

ti = -2:.05:2; 
[XI,YI] = meshgrid(ti,ti);
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