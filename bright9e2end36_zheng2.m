%% mingchang 9 xuan 2 end 72 zheng
%% caculate all 9 led then next different round  9to8
%%
%clc;
%clear;
%close all;

loop=1;
alpha = eps;%1
beta = 1e8;%1e3
gamma = 8;

%% ����
addpath(['./FP_Func']);
addpath('natsortfiles');


LEDheight=67.5;%90   67.5
%dd = sqrt((-hhled*ds_led-img_center(1)).^2+(-vvled*ds_led-img_center(2)).^2+z_led.^2);


LEDgap = 4;%LED���Ϊ4mm
waveLength = 0.6292;%%������Ĳ���
m = 384; 
n = 384;
%  �߷ֱ��������ͼ��ߴ� [m,n]�ֱ�������������������
m1 = m/3; %
n1 = n/3; %���������ͼ���С             n = 3*n1


mag = 8.1485;
dpix_c = 6.5; %6.5um pixel size on the sensor plane
dpix_m = dpix_c/mag; 
dkx = 1/(m1*dpix_m);
dky = 1/(m1*dpix_m);
NA = 0.1;
% maximum spatial frequency set by NA
um_m = NA/waveLength;
rx = 1:m1;
[rxx,ryy] = meshgrid(rx-round((m1+1)/2));
ridx = sqrt(rxx.^2+ryy.^2);
um_idx = um_m/dkx;
% assume a circular pupil function, lpf due to finite NA
CTF = double(ridx<um_idx);

nstart = [973,1173];
%nstart = [800,1400];
%nstart = [1017,1217];
ncent = [1080,1280];
%img_ncent = nstart-ncent+m1/2;
img_center = (nstart-ncent+m1/2)*dpix_m;
xlo = struct2array(load('xloylo.mat', 'xlo'));%293��ˮƽ����
ylo = struct2array(load('xloylo.mat', 'ylo'));%293����ֱ����

dd2 = sqrt((-xlo-img_center(1)).^2+(-ylo-img_center(2)).^2+LEDheight.^2);
% kxmid = round((-xlo-img_center(1))./dd2/lambda/du);
% kymid = round((-ylo-img_center(2))./dd2/lambda/du);

% kx_relative = (xlo+img_center(2)/1000)./dd2;
% ky_relative = (ylo+img_center(1)/1000)./dd2;

kx_relative = sin(atan((xlo-img_center(1)/1000)/LEDheight));
ky_relative = sin(atan((ylo-img_center(2)/1000)/LEDheight));
kx_relative = reshape(kx_relative,1,293);
ky_relative = reshape(ky_relative,1,293);

kx=kx_relative/waveLength;
ky=ky_relative/waveLength;
k1=kx;
k2=ky;
ky=k1;
kx=k2;

%newimSeqLowRes=zeros(m1,m1,37);%Ԥ�����ڴ�
ledidx = struct2array(load('expt_lit-8.mat', 'ledidx')) ;
%ledidx = ledidx(1:592);
%ledidx(1,173) = 203;
%usedled = sort(ledidx);
sort9 = reshape(ledidx,8,[]);
sort9T = sort9.';
sort9T1 = sort9T(1:37,:);
sort9T2 = sort9T(37:74,:);
sort9T3 = sort9T(74:110,:);
sort9T4 = sort9T(110:147,:);
sort9T5 = sort9T(147:184,:);
sort9T6 = sort9T(184:220,:);
sort9T7 = sort9T(220:257,:);
sort9T8 = sort9T(257:293,:);

%% read in all images into the memory first
filedir1 = '/home/jnu733/qwm/data/usaf1led/';
%filedir1 = '/home/jnu733/qwm/data/cell8led/';
%filedir2 = '/home/jnu733/qwm/data/cell1led/';

% D:/workapp/matlab/workstation/QWM-FPM/data/8led/tif/
% /home/jnu733/qwm/data
imglist1 = dir([filedir1,'*.tif']);
N = natsortfiles({imglist1.name});
%numlit = 8;
%nn1 = 128; nn2 = 128;

fprintf('loading the images...\n');
tic;
Nimg = length(imglist1);
Iall = zeros(2160,2560,Nimg);

Ibk = zeros(Nimg,1);
for mm1 = 1:Nimg
    fn = [filedir1,N{mm1}];  %added for sorting purpose
    disp(fn);
    I = double(imread(fn));

    bk1 = mean2(double(I(1:100,1:100)));
    bk2= mean2(double(I(492:600,2380:2520)));
    Ibk(mm1) = mean([bk1,bk2]);
    
    
    Iall(:,:,mm1) = I  ;

       
%      if Ibk(mm1)>1000 && mm1 >2
%             Ibk(mm1) = Ibk(mm1-1);
%      end
     if Ibk(mm1)>1000
            Ibk(mm1) = Ibk(mm1-1);
     end
     
     
end
% Ibk(1) = Ibk(3);
% Ibk(2) = Ibk(3);

for mm1 = 1:Nimg
    Itmp = Iall(:,:,mm1);
    Itmp = Itmp-Ibk(mm1);
%     Itmp = awgn(Itmp,0,'measured');
    Itmp(Itmp<0) = 0;
    Iall(:,:,mm1) = Itmp;
    
end

Iallend = double(Iall(nstart(1):nstart(1)+m1-1,nstart(2):nstart(2)+m1-1,:));
clear I Iall


newimSeqLowRes = double(Iallend);


% num_list = setdiff(1:293,[127,128,129,146,147,148,165,166,167]);
% usedled = sort(num_list);
% num_list = reshape(num_list(randperm(length(num_list))),4,[])';
% newimSeqLowRes2 = zeros(128,128,71);
% for k = 1:71
%     newimSeqLowRes2(:,:,k) = newimSeqLowRes(:,:,num_list(k,1))+newimSeqLowRes(:,:,num_list(k,2))...
%                             +newimSeqLowRes(:,:,num_list(k,3))+newimSeqLowRes(:,:,num_list(k,4));
% end    

bright = [127,128,129,146,147,148,165,166,167];
% num_list = nchoosek(bright,2);
num_list2 = struct2array(load("num_list2.mat","num_list2"));
newimSeqLowRes2 = zeros(128,128,36);
for k = 1:36
    newimSeqLowRes2(:,:,k) = newimSeqLowRes(:,:,num_list2(k,1))+newimSeqLowRes(:,:,num_list2(k,2));
end  
%save("num_list2.mat","num_list2")
%save("seq288.mat","num_list","-append")

%objectlow = newimSeqLowRes(:,:,1);
%I_one =(newimSeqLowRes(:,:,147));
objectlow = sqrt(abs(newimSeqLowRes2(:,:,1)));
figure(86);
subplot(221); imagesc(abs(objectlow)); axis image; colormap gray; colorbar;
title('ampl(o)');
subplot(222); imagesc(angle(objectlow)); axis image; colormap gray; colorbar;
title('phase(o)');
subplot(223); imagesc(angle(objectlow)); axis image; colormap gray; colorbar;
title('ampl(O)');
subplot(224); imagesc(angle(objectlow)); axis image; colorbar;
title('phase(O)');




%% �����ǻָ�����
%seq = gseq(arraysize);%����ָ���˳�򣬴�����()����Ե
seq293 = struct2array(load('seq288.mat', 'seq293'));%��������˳��
% seq288(45) = [];
% seq288 = [seq288;293];
% objectRecover = imresize(sqrt(newimSeqLowRes(:,:,1)),3,'bilinear');%�Ը߷ֱ�������ĳ�ʼ�²�(ȫ1����)
% objectRecoverFT = fftshift(fft2(objectRecover));%�ѳ�ʼ�²�ת������Ҷ��(�ڸ���Ҷ�򴴽���һ������)
F = @(x) fftshift(fft2(x));
Ft = @(x) ifft2(ifftshift(x));
upsamp = @(x) padarray(x,[(m-m1)/2,(m-m1)/2]);
objectRecoverFT = F(sqrt(newimSeqLowRes2(:,:,1)));
objectRecoverFT = upsamp(objectRecoverFT);

objectRecoverFT0 = objectRecoverFT;

pupil0=double(CTF);
pupil=pupil0;





fprintf('| iter |  rmse    |\n');
for j=1:20, fprintf('-'); end
fprintf('\n');
err2 = 0;
err = [];
iter = 0;

fprintf('| %2d   | %.2e |\n',iter,err2);
for tt=1:loop
    

    

    %% nine led zheng1
    for i4 = 1:8
            iter = iter+1;
            err2 = 0;
            gamma = gamma/2;
            for i3=1:9 
                i2=seq293(i3);

                [x,y]=find(num_list2==i2);


                    lowResFT = zeros(m1,m1,2);
                    im_lowRes = zeros(m1,m1,2);
                    for r=1:2
                        kxc=round((n+1)/2+kx(1,num_list2(x(i4),r))/dkx);%round�Ƕ������ڵ�Ԫ��ȡ��
                        kyc=round((m+1)/2+ky(1,num_list2(x(i4),r))/dky);
                        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
                        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
                        lowResFT(:,:,r)=(m1/m)^2*objectRecoverFT(kyl:kyh,kxl:kxh).*pupil;
                        im_lowRes(:,:,r) = abs(ifft2(ifftshift(lowResFT(:,:,r))));
                    end

                    %im_lowRes_pinpu=lowResFT_1+lowResFT_2+lowResFT_3;
                    im_lowRes_oI=(abs(im_lowRes(:,:,1))).^2+(abs(im_lowRes(:,:,2))).^2;

                    kxc=round((n+1)/2+kx(1,i2)/dkx);%round�Ƕ������ڵ�Ԫ��ȡ��
                    kyc=round((m+1)/2+ky(1,i2)/dky);
                    kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
                    kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
                    lowResFT_ori=(m1/m)^2*objectRecoverFT(kyl:kyh,kxl:kxh).*pupil;
                    im_lowRes_ori = ifft2(ifftshift(lowResFT_ori));

                    %im_lowRes = (m/m1)^2*imSeqLowRes(:,:,i2).*exp(1i.*angle(im_lowRes));%�滻���
                    %im_lowRes = (m/m1)^2*sqrt(abs(newimSeqLowRes(:,:,x))./max(max(abs(newimSeqLowRes(:,:,x)))))./sqrt(abs(im_lowRes_oI)./max(max(abs(im_lowRes_oI)))).*im_lowRes_ori;%�滻���
                    im_lowResk = (m/m1)^2*sqrt(abs((newimSeqLowRes2(:,:,i4))))./sqrt(abs(im_lowRes_oI)).*im_lowRes_ori;%�滻���

                    lowResFT=fftshift(fft2(im_lowResk)).*pupil;%ת��Ϊ����Ҷ��
                    OP_diff=(m1/m)^2*lowResFT-lowResFT_ori;
                    OFTdown = objectRecoverFT(kyl:kyh,kxl:kxh);
                    %objectRecoverFT(kyl:kyh,kxl:kxh)=(1-pupil).*objectRecoverFT(kyl:kyh,kxl:kxh)+0.1*lowResFT;%�ٷ���ȥ���³�ʼ�²�ĸ���Ҷ������
                    %objectRecoverFT(kyl:kyh,kxl:kxh)=objectRecoverFT(kyl:kyh,kxl:kxh)+1*abs(pupil).*(lowResFT-lowResFT_ori./(pupil+eps));
                    objectRecoverFT(kyl:kyh,kxl:kxh)=objectRecoverFT(kyl:kyh,kxl:kxh)+gamma*abs(pupil).*conj(pupil).*(OP_diff)./max(max(abs(pupil)))./(abs(pupil).^2+alpha);
                    pupil =  pupil+gamma/max(max(abs(objectRecoverFT))).*(abs(OFTdown).*conj(OFTdown)).*(OP_diff)./(abs(OFTdown).^2+beta).*pupil0;

                    err2 = err2+sqrt(sum(sum(((abs(newimSeqLowRes2(:,:,i4))-abs(im_lowRes_oI)).^2))));
             end

     



    fprintf('| %2d   | %.2e |\n',iter,err2);
    objectRecover=ifft2(ifftshift(objectRecoverFT));
    figure(38);
    subplot(221); imagesc(abs(objectRecover)); axis image; colormap gray; colorbar;
    title('ampl(o)');
    subplot(222); imagesc(angle(objectRecover)); axis image; colormap gray; colorbar;
    title('phase(o)');
    subplot(223); imagesc(abs(log((objectRecoverFT+1)))); axis image; colormap gray; colorbar;
    title('objectRecoverFT');
    subplot(224); imagesc(abs(log(objectRecoverFT0+1))); axis image; colormap gray; colorbar;
    title('objectRecoverFT0');%pupil    log(objectRecoverFT0)
    end

end
objectRecover=ifft2(ifftshift(objectRecoverFT));
figure(37);
subplot(221); imagesc(abs(objectRecover)); axis image; colormap gray; colorbar;
title('ampl(o)');
subplot(222); imagesc(angle(objectRecover)); axis image; colormap gray; colorbar;
title('phase(o)');
subplot(223); imagesc(abs((pupil))); axis image; colormap gray; colorbar;
title('abs((pupil)');
subplot(224); imagesc(abs(angle((pupil)))); axis image; colormap gray; colorbar;
title('angle((pupil))');%pupil    log(objectRecoverFT0)
% (log(objectRecoverFT)),[]);title('objectRecoverFT');
% (log(objectRecoverFT0)),[]);title('objectRecoverFT0');
% [ms,ps]=mse_psnr(objAmplitude,abs(objectRecover),8);
fprintf('processing completes\n');

%ampl=abs(objectRecover);
%am_high = max(max(ampl));
%am_low = min(min(ampl));
%ampl = (ampl-am_low)./(am_high-am_low);
%imwrite(ampl,'result.png');
%savePath = ['D:/workapp/matlab/workstation/QWM-FPM/data/8ledsimulate/result/'];
ampl=abs(objectRecover);
am_high = max(max(ampl));
am_low = min(min(ampl));
ampl = (ampl-am_low)./(am_high-am_low);
imwrite(ampl,'/home/jnu733/qwm/result/matlab_result/2led/usaf2ampli_36_zheng9to8.png','png');
% fileName = ['simulate293ampli.tif'];
% fullPath = fullfile(savePath,fileName); % ��ȡ�����ı���·�����ļ���
% imwrite(I_est,fullPath, 'tif', 'Compression', 'none'); % ����Ϊ16λ��ͨ��tifͼƬ���ر�ѹ��

ampl=angle(objectRecover);
am_high = max(max(ampl));
am_low = min(min(ampl));
ampl = (ampl-am_low)./(am_high-am_low);
imwrite(ampl,'/home/jnu733/qwm/result/matlab_result/2led/usaf2angle_36_zheng9to8.png','png');
% fileName = ['simulate293angle.tif'];
% fullPath = fullfile(savePath,fileName); % ��ȡ�����ı���·�����ļ���
% imwrite(I_est,fullPath, 'tif', 'Compression', 'none'); % ����Ϊ16λ��ͨ��tifͼƬ���ر�ѹ��



