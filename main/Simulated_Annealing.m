% =========================================================================
% Introduction
% =========================================================================
% 模拟退火算法不适用于找高维成本函数最小值对应的全局最优解，考虑解的变化区
% 间都在（0，1），则1维变量x中每个温度需要搜索足够的点才能有效覆盖解（这取
% 决于变量的精度位数），持续降温后可以得到全局最优解。但是对于25*25=625维实
% 数空间，假如每个维度搜索100个点，需要搜索的总点数就是100^625个，因此计算
% 量是非常巨大的。

% Author: Ya-lowkey (cldeng881@gmail.com)
% =========================================================================
clc
clear
close
addpath(genpath('./imgs'))
addpath(genpath('./function'))
% pic1=imread('cameraman.bmp');
  pic1=ones(100,100);
  pic1(25:75,25:75)=zeros(51,51);
pic1=double(pic1);
[row1,col1]=size(pic1);
imshow(pic1,[])
%sz=100;%设置模拟退火算法搜索图片的变量维度大小
%pic1=double(pic1(floor(row1/2)-sz/2:floor(row1/2)+sz/2-1,floor(col1/2)-sz/2:floor(col1/2)+sz/2-1));%调节剪裁大小
[row,col]=size(pic1);
pic1=pic1;
wavelen=532e-9;
dist=1;
pixsize=4e-6;

T0 = 100;   % 初始温度
T = T0; % 迭代初始温度T0
epoch =350;  % 最大迭代次数

Lk = 100;  % 每个温度下的迭代次数
alfa =0.95;  % 温度衰减系数

% 随机生成一个初始解
x0=rand(row,col);  
pic0=exp(1i*2*pi*x0);
I=propagate(pic0,dist,pixsize,wavelen);
y0=sum((I.*conj(I)-pic1).^2,'all')./(row*col);
min_y = y0;     
% MINY = zeros(epoch,1); 
% 模拟退火过程
%%
for iter = 1 : epoch  % 总迭代次数
    iter
    for i = 1 : Lk  % 每个温度下迭代Lk次
        
        z = (rand(row,col)-rand(row,col)); % 均匀分布（-1，1）
        x_new = x0+ z*(T./T0); % 随机更新产生不超过[-1,1]区间的新坐标
        x_new(x_new>1)=1;
        x_new(x_new<0)=0;
        
        x1 = x_new;    
           pic0=exp(1i*2*pi*x1./2);
           I=propagate(pic0,dist,pixsize,wavelen);
        y1=sum((I.*conj(I)-pic1).^2,'all')./(row*col);% 计算均方差
        
        if y1 < y0    % 误差降低自然接受新点
            x0 = x1; 
            y0 = y1;
%         else
%             p = exp(-(y1-y0)/T); % 误差升高按概率p接受新坐标点
%             if rand(1) < p   % 此条件满足则接受新坐标点，条件会随着温度降低变得很难满足，导致坐标维持原坐标点
%                 x0 = x1; % 更新当前解为新解
%                 y0 = y1;
%             end
        end
        
        if y0 < min_y  % 通过前后比较，存储该温度下的最小值
            min_y = y0;  
            best_x = x0; 
        end
    end
    if iter>50
        alfa=0.99; % 温度按衰减率下降，逐渐接近0
        Lk=150;
    end  
     if iter>300
         Lk=350;
         alfa=0.999;
     end
    MINY(iter) = min_y; % 保存每个温度迭代下得到的最小值误差
    T=T*alfa
end
%%
subplot(2,1,1)
plot(1:iter,MINY,'b-')
    xlabel('iteration')
    ylabel('MSE')
    axis([0,epoch,0,1.2])
subplot(2,1,2)
imshow(abs(I).^2,[])
title('衍射效果')
set(gcf,'color','w')















































