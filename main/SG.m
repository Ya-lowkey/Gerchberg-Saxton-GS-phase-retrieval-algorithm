%%
%--------------------------------------------------------------------------
% Author: Ya-lowkey (cldeng881@gmail.com)
% ��ϸ�Ƶ����ע΢�Ź��ں� @���ӿ���
%--------------------------------------------------------------------------
%%
%�ݶ��½�ʵ����λ�ָ�
clc
clear
close all
addpath(genpath('./imgs'))
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
for j=1:100
A=fft2(aphi);
ab=A.*conj(A);
ab=ab./max(ab(:));
dA=ab-I2;
%dphi=real(propagate(dA.*A.*2./(r*c),-dist,pixsize,wavelen).*(-1i*conj(aphi)));%��������ʱѡ������ݶȻش�
dphi=real(ifft2(dA.*A.*2./(r*c)).*(-1i*conj(aphi)));
g=dphi;

%��������Ӧadam�Ż�
v=beta1*v+(1-beta1)*g;
s=beta2*s+(1-beta2)*g.^2;
vc=v./(1-beta1^j);
sc=s./(1-beta2^j);
m=vc./(sqrt(sc)+b).*r;
phi=phi-m;

aphi=exp(1i*phi);
subplot(1,2,1)
imshow(abs(fft2(aphi)).^2,[])
title('ÿ�ε���������Ч��')
subplot(1,2,2)
plot(j,sum(dA(:).^2)./(r*c),'b.')
xlabel('iteration')
ylabel('MSE')
hold on
drawnow
end
%imshow(I2)

