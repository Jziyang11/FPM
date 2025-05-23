%% 模拟前向模型过程
%% 生成高分辨率物体作为输入
objectAmplitude = double(imread('cameraman.tif'));
phase = double(imread('westconcordorthophoto.png'));
phase = pi*imresize(phase,[256 256])./max(max(phase));   % 归一化相位到[0,pi]
object = objectAmplitude.*exp(1i.* phase);    % 复振幅图像
subplot(2,3,1)
imshow(abs(object),[]);title('输入高分辨率物体振幅图像')
% imshow用于显示图像，输入为实数矩阵，对于复数矩阵默认输出幅度。
% 且用于模拟传感器在成像过程中丢失相位信息

%% 定义光学系统参数
waveLength = 0.63e-6;
k0 = 2*pi/waveLength;
spsize = 2.75e-6;     % ccd成像像素尺寸 
psize = spsize / 4;   % 重建图像像素尺寸
NA = 0.08;

%% 生成LED灯照明波矢
arraysize = 15;
xlocation = zeros(1,arraysize^2);  % 坐标x：长度为15^2的一维行向量
ylocation = zeros(1,arraysize^2);
LEDgap = 4;     % 相邻led灯间距
LEDheight = 90; % led到样本的距离

for i=1:arraysize   % 从左上到右下
    xlocation(1,1+arraysize*(i-1):15+arraysize*(i-1)) = (-(arraysize-1)/2:1:(arraysize-1)/2)*LEDgap;
    ylocation(1,1+arraysize*(i-1):15+arraysize*(i-1)) = ((arraysize-1)/2-(i-1))*LEDgap;
end
kx_relative = -sin(atan(xlocation/LEDheight));  % x方向波矢
ky_relative = -sin(atan(ylocation/LEDheight));  % y方向波矢

%% 模拟低分辨率图像采集
[m,n] = size(object); % 输入高分辨率图像尺寸
m1 = m/(spsize/psize);n1 = n/(spsize/psize); % 输出低分辨率图像尺寸（像素数）
imSeqLowRes = zeros(m1, n1, arraysize^2); % 存储输出低分辨率图像序列 225个m1×n1矩阵
kx = k0 * kx_relative;
ky = k0 * ky_relative;  % 绝对波矢
dkx = 2*pi/(psize*n);   % 频率采样间隔(x方向)
dky = 2*pi/(psize*m);  
cutoffFrequency = NA * k0;
kmax = pi/spsize;
[kxm,kym] = meshgrid(-kmax:kmax/((n1-1)/2):kmax,-kmax:kmax/((n1-1)/2):kmax);
CTF = ((kxm.^2+kym.^2)<cutoffFrequency^2);  % 无像差的光瞳函数
% 有像差，定义有离焦像差的光瞳函数，离焦距离z
z = 10e-6; kzm = sqrt(k0^2-kxm.^2-kym.^2);
pupil = exp(1i.*z.*real(kzm)).*exp(-abs(z).*abs(imag(kzm)));
aberratedCTF = pupil.*CTF;  

%close all;
clc; % 清除命令行
objectFT = fftshift(fft2(object));
for tt =1:arraysize^2
    kxc = round((n+1)/2+kx(1,tt)/dkx); 
    kyc = round((m+1)/2+ky(1,tt)/dky); % 计算频域中心位置
    kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
    kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2); % 确定频域截取范围
    imSeqLowFT = (m1/m)^2 * objectFT(kyl:kyh,kxl:kxh).*aberratedCTF;  % 经过低通滤波
    imSeqLowRes(:,:,tt) = abs(ifft2(ifftshift(imSeqLowFT))); % 只取强度信息
end
subplot(2,3,2)
imshow(imSeqLowRes(:,:,1),[]);title('低分辨率图像')
% close all;

%% 恢复高分辨率图像
seq = gseq(arraysize);   % 定义恢复的顺序
objectRecover = ones(m,n); % 初始化猜测
objectRecoverFT = fftshift(fft2(objectRecover));
loop = 5; % 定义迭代次数
for tt=1:loop
    for i3=1:arraysize^2
        i2=seq(i3);
        kxc = round((n+1)/2+kx(1,i2)/dkx);
        kyc = round((m+1)/2+ky(1,i2)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        lowResFT = (m1/m)^2 * objectRecoverFT(kyl:kyh,kxl:kxh).*aberratedCTF;
        im_lowRes = ifft2(ifftshift(lowResFT));
        im_lowRes = (m/m1)^2 * imSeqLowRes(:,:,i2).*exp(1i.*angle(im_lowRes)); 
        lowResFT=fftshift(fft2(im_lowRes)).*CTF.*(1./pupil); % 反转光瞳函数，补偿像差
        objectRecoverFT(kyl:kyh,kxl:kxh)=(1-CTF).*objectRecoverFT(kyl:kyh,kxl:kxh) + lowResFT;                   
    end
end
objectRecover=ifft2(ifftshift(objectRecoverFT));
subplot(2,3,3)
imshow(abs(objectRecover),[]);title('恢复的幅度')
subplot(2,3,4)
imshow(angle(objectRecover),[]);title('恢复的相位')
subplot(2,3,5)
imshow(log(objectRecoverFT),[]);title('恢复的频谱')
