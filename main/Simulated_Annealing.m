% =========================================================================
% Introduction
% =========================================================================
% ģ���˻��㷨���������Ҹ�ά�ɱ�������Сֵ��Ӧ��ȫ�����Ž⣬���ǽ�ı仯��
% �䶼�ڣ�0��1������1ά����x��ÿ���¶���Ҫ�����㹻�ĵ������Ч���ǽ⣨��ȡ
% ���ڱ����ľ���λ�������������º���Եõ�ȫ�����Ž⡣���Ƕ���25*25=625άʵ
% ���ռ䣬����ÿ��ά������100���㣬��Ҫ�������ܵ�������100^625��������ĳ�
% �����ģ����һ������ͼƬ��ά�ȼ����ؿ������е��ڡ���˼������Ƿǳ��޴�ġ�

% Author: Ya-lowkey (cldeng881@gmail.com)
% =========================================================================
clc
clear
addpath(genpath('./imgs'))
addpath(genpath('./function'))
pic1=imread('boat.bmp');
pic1=double(pic1(:,:,1));
[row1,col1]=size(pic1);
sz=25;%����ģ���˻��㷨����ͼƬ�ı���ά�ȴ�С���������25*25=625��ά��
pic1=double(pic1(floor(row1/2)-sz/2:floor(row1/2)+sz/2-1,floor(col1/2)-sz/2:floor(col1/2)+sz/2-1));%���ڼ��ô�С
[row,col]=size(pic1);
imshow(pic1,[])
%%
pic1=pic1./250;
wavelen=532e-9;
dist=1;
pixsize=4e-6;

T0 = 100;   % ��ʼ�¶�
T = T0; % ������ʼ�¶�T0
epoch = 5;  % ����������
Lk = 10;  % ÿ���¶��µĵ�������
alfa = 0.95;  % �¶�˥��ϵ��

% �������һ����ʼ��
x0=rand(row,col);  
pic0=exp(1i*2*pi*x0);
I=propagate(pic0,dist,pixsize,wavelen);
I=abs(I);
ae=corrcoef(I,pic1)
y0 = ae(1,2);

min_y = y0;     
MINY = zeros(epoch,1); 

% ģ���˻����
for iter = 1 : epoch  % �ܵ�������
    iter
    for i = 1 : Lk  % ÿ���¶��µ���Lk��
        
        z = rand(row,col)-rand(row,col); % ��̫�ֲ��ı���
        x_new = x0+ z.*(T./T0); % ������²���������[-1,1]�����������
       [h1,~]=size(x_new(x_new<0));
       [h2,~]=size(x_new(x_new>1));
       x_new(x_new<0)=x_new(x_new<0)+ones(h1,1);%��������Χ��������������ۻ���������1.5���0.5��-0.5���0.5
       x_new(x_new>1)=x_new(x_new>1)-ones(h2,1);
       
        x1 = x_new;    
           pic0=exp(1i*2*pi*x1);
           I=propagate(pic0,dist,pixsize,wavelen);
           ae=corrcoef(I,pic1);% �����������ص�������ֵ�͵�ƽ��ֵ
        y1 = ae(1,2);  
        if y1 > y0    % ������Ȼ�����µ�
            x0 = x1; 
            y0 = y1;
        else
            p = exp(-(y1- y0)/T); % ������߰�����p�����������
            if rand(1) < p   % ���������������������㣬�����������¶Ƚ��ͱ�ú������㣬��������ά��ԭ�����
                x0 = x1; % ���µ�ǰ��Ϊ�½�
                y0 = y1;
            end
        end
        
        if y0 > min_y  % ͨ��ǰ��Ƚϣ��洢���¶��µ���Сֵ
            min_y = y0;  
            best_x = x0;  
        end
    end
    MINY(iter) = min_y; % ����ÿ���¶ȵ����µõ�����Сֵ���
    T = alfa*T;   % �¶Ȱ�˥�����½����𽥽ӽ�0
   
end
%�������ϵ�����ߣ�Խ�ӽ�1����ṹ���ƶ�Խ��
 plot(1:epoch,MINY,'b*-')
 xlabel('iteration')
 ylabel('ae')

 %% 
 pic0=exp(1i*2*pi*best_x);
 I=propagate(pic0,dist,pixsize,wavelen);
 I=abs(I)./max(max(abs(I)));
 I=imresize(I,[1024,1024]);
 imshow(I)















































