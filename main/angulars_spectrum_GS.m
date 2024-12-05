
% Author: Ya-lowkey (cldeng881@gmail.com)
% =========================================================================
clc
clear 
addpath(genpath('./imgs'))
addpath(genpath('./function'))
I1=double(imread('boat.bmp','bmp'));%ÊäÈëÍ¼Ïñ1£¨²Î¿¼Í¼Ïñ£©
I1=I1(:,:,1);
[row,col]=size(I1);
wavelen=532e-9;
dist=1;
pixsize=4e-6;
figure(1)
imshow(I1,[]);%
figure(2)
meshc(I1)%nint8 (0,255)
%%
epoch=50;
aphi=exp(1i*2*pi.*rand(row,col));% Forward diffraction
for i=1:epoch
    
I=propagate(aphi,dist,pixsize,wavelen);
R_cor=corrcoef(I1,abs(I));

R1=fft2(I1);
R1=R1./max(max(R1));
Inorm=255.*abs(I)./max(max(abs(I)));
R2=fft2(Inorm);
R2=R2./max(max(R2));
R=ifftshift(ifft2(R1.*conj(R2)));% Cross-correlation operations are used to evaluate the structural similarity of images
R=R./max(max(R));

Iabs=abs(I)./max(max(abs(I)));
phi0=(angle(I)+pi)./(2*pi);% Obtain the phase angle of the computational graph I and map it to the interval (0, 1).
aphi0=exp(1i*2*pi.*phi0);
I_target=(I1+rand(1,1).*(I1-Iabs)).*aphi0;% Amplitude constraint of the target object
I_start=propagate(I_target,-dist,pixsize,wavelen);% Back diffraction
phi=(angle(I_start)+pi)./(2*pi);
aphi=exp(1i*2*pi.*phi);

subplot(2,2,1);
imshow(abs(I),[])
title('Diffraction pattern')
subplot(2,2,2);
imshow(abs(phi))%imshow phase
title('phase of DOE')
subplot(2,2,3);
imshow(abs(R),[])
title('Cross-correlation')
subplot(2,2,4);
plot(i,R_cor(1,2),'b.')
hold on
axis([0,epoch,0,2])
title('corrcoef')

pause(0.01)
end
%%
II=abs(I)./max(max(abs(I)));
imshow(uint8(II.*255))
imwrite(uint8(II.*255),'result_agul.bmp','bmp');








