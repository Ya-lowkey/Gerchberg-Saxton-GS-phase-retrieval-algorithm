%%
%--------------------------------------------------------------------------
% Author: Ya-lowkey (cldeng881@gmail.com)
% 详细推导请关注微信公众号 @智子科普
%--------------------------------------------------------------------------
%%
%梯度下降实现相位恢复
clc
clear
close all
addpath(genpath('./imgs'))
addpath(genpath('./function'))
I2=double(imread('lake.bmp','bmp'));
I2=I2./max(I2(:));
[r,c]=size(I2);
phi=2*pi*rand(r,c);
aphi=exp(1i*phi);
beta1=0.9;
beta2=0.999;
v=zeros(r,c);
s=zeros(r,c);
b=0.0001;
r=1;

wavelen=532e-9;%波长m
dist=1;%衍射距离m
pixsize=4e-6;%像素尺寸m

for j=1:100
A=fft2(aphi);
%A=propagate(aphi,dist,pixsize,wavelen);%仿真衍射时选择这个正向衍射
ab=A.*conj(A);
ab=ab./max(ab(:));
dA=ab-I2;
%dphi=real(propagate(dA.*A.*2./(r*c),-dist,pixsize,wavelen).*(-1i*conj(aphi)));%仿真衍射时选择这个梯度回传
dphi=real(ifft2(dA.*A.*2./(r*c)).*(-1i*conj(aphi)));
g=dphi;

%启用自适应adam优化
v=beta1*v+(1-beta1)*g;
s=beta2*s+(1-beta2)*g.^2;
vc=v./(1-beta1^j);
sc=s./(1-beta2^j);
m=vc./(sqrt(sc)+b).*r;
phi=phi-m;

aphi=exp(1i*phi);
subplot(1,2,1)
imshow(abs(A).^2,[])
title('每次迭代的衍射效果')
subplot(1,2,2)
plot(j,sum(dA(:).^2)./(r*c),'b.')
xlabel('iteration')
ylabel('MSE')
hold on
drawnow
end
%imshow(I2)
