clc; clear all; close all;

x=-8:0.5:8;

[XX,YY]=meshgrid(x);

r=sqrt(XX.^2+YY.^2)+eps;

Z=sin(r)./r;

surf(Z);

theAxes=axis;

fmat=moviein(20);

for j=1:20

surf(sin(2*pi*j/20)*Z,Z)

axis(theAxes)

fmat(:,j)=getframe;

end

movie(fmat,10)