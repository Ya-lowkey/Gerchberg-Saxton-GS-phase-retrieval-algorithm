% =========================================================================
% Introduction
% =========================================================================
% 模拟退火算法不适用于找高维成本函数最小值对应的全局最优解，考虑解的变化区
% 间都在（0，1），则1维变量x中每个温度需要搜索足够的点才能有效覆盖解（这取
% 决于变量的精度位数），持续降温后可以得到全局最优解。但是对于25*25=625维实
% 数空间，假如每个维度搜索100个点，需要搜索的总点数就是100^625个，下面的程
% 序就是模拟这一场景，图片的维度即像素可以自行调节。因此计算量是非常巨大的。

% Author: Ya-lowkey (cldeng881@gmail.com)
% =========================================================================
clc
clear
addpath(genpath('./imgs'))
addpath(genpath('./function'))
pic1=imread('boat.bmp');
pic1=double(pic1(:,:,1));
[row1,col1]=size(pic1);
sz=25;%设置模拟退火算法搜索图片的变量维度大小，这里就是25*25=625个维度
pic1=double(pic1(floor(row1/2)-sz/2:floor(row1/2)+sz/2-1,floor(col1/2)-sz/2:floor(col1/2)+sz/2-1));%调节剪裁大小
[row,col]=size(pic1);
imshow(pic1,[])
%%
pic1=pic1./250;
wavelen=532e-9;
dist=1;
pixsize=4e-6;

T0 = 100;   % 初始温度
T = T0; % 迭代初始温度T0
epoch = 5;  % 最大迭代次数
Lk = 10;  % 每个温度下的迭代次数
alfa = 0.95;  % 温度衰减系数

% 随机生成一个初始解
x0=rand(row,col);  
pic0=exp(1i*2*pi*x0);
I=propagate(pic0,dist,pixsize,wavelen);
I=abs(I);
ae=corrcoef(I,pic1)
y0 = ae(1,2);

min_y = y0;     
MINY = zeros(epoch,1); 

% 模拟退火过程
for iter = 1 : epoch  % 总迭代次数
    iter
    for i = 1 : Lk  % 每个温度下迭代Lk次
        
        z = rand(row,col)-rand(row,col); % 正太分布的倍率
        x_new = x0+ z.*(T./T0); % 随机更新产生不超过[-1,1]区间的新坐标
       [h1,~]=size(x_new(x_new<0));
       [h2,~]=size(x_new(x_new>1));
       x_new(x_new<0)=x_new(x_new<0)+ones(h1,1);%将超出范围的坐标点周期性折回来，比如1.5变成0.5，-0.5变成0.5
       x_new(x_new>1)=x_new(x_new>1)-ones(h2,1);
       
        x1 = x_new;    
           pic0=exp(1i*2*pi*x1);
           I=propagate(pic0,dist,pixsize,wavelen);
           ae=corrcoef(I,pic1);% 计算所有像素点误差绝对值和的平均值
        y1 = ae(1,2);  
        if y1 > y0    % 误差降低自然接受新点
            x0 = x1; 
            y0 = y1;
        else
            p = exp(-(y1- y0)/T); % 误差升高按概率p接受新坐标点
            if rand(1) < p   % 此条件满足则接受新坐标点，条件会随着温度降低变得很难满足，导致坐标维持原坐标点
                x0 = x1; % 更新当前解为新解
                y0 = y1;
            end
        end
        
        if y0 > min_y  % 通过前后比较，存储该温度下的最小值
            min_y = y0;  
            best_x = x0;  
        end
    end
    MINY(iter) = min_y; % 保存每个温度迭代下得到的最小值误差
    T = alfa*T;   % 温度按衰减率下降，逐渐接近0
   
end
%绘制相关系数曲线，越接近1代表结构相似度越高
 plot(1:epoch,MINY,'b*-')
 xlabel('iteration')
 ylabel('ae')

 %% 
 pic0=exp(1i*2*pi*best_x);
 I=propagate(pic0,dist,pixsize,wavelen);
 I=abs(I)./max(max(abs(I)));
 I=imresize(I,[1024,1024]);
 imshow(I)















































