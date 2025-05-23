%% tianlei 1led liang && 4led anchang
%clc;
%clc;
clear all;
%close all;

loop=5;
alpha = eps;%1
beta = eps;%1e16;%1e20;%1e3
%gamma1 =8;
gamma2 =0.002;%object
gamma3=0.05;%pupil
num_image = 45;

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

ky=kx_relative/waveLength;
kx=ky_relative/waveLength;

%% read in all images into the memory first
filedir1 = '/home/jnu733/qwm/data/cell1led/1led/';
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

       
 %     if Ibk(mm1)>6000 && mm1 >2
 %            Ibk(mm1) = Ibk(mm1-1);
  %    end
%      if Ibk(mm1)>1000&& mm1 >2
%             Ibk(mm1) = Ibk(mm1-1);
%      end
%      
     
end
% Ibk(1) = Ibk(3);
% Ibk(2) = Ibk(3);
% Ibk(:)= 8000;
for mm1 = 1:Nimg
    Itmp = Iall(:,:,mm1);
    Itmp = Itmp-Ibk(mm1);
%     Itmp = awgn(Itmp,0,'measured');
   % Itmp(Itmp<0) = 0;
    Iall(:,:,mm1) = Itmp;
    
end

Iallend = double(Iall(nstart(1):nstart(1)+m1-1,nstart(2):nstart(2)+m1-1,:));
%clear I Iall


newimSeqLowRes = double(Iallend);


% num_list = setdiff(1:293,[127,128,129,146,147,148,165,166,167]);
% usedled = sort(num_list);
% num_list = reshape(num_list(randperm(length(num_list))),4,[])';
% newimSeqLowRes2 = zeros(128,128,71);
% for k = 1:71
%     newimSeqLowRes2(:,:,k) = newimSeqLowRes(:,:,num_list(k,1))+newimSeqLowRes(:,:,num_list(k,2))...
%                             +newimSeqLowRes(:,:,num_list(k,3))+newimSeqLowRes(:,:,num_list(k,4));
% end    


% num_list = nchoosek(bright,2);
% num_list2 = struct2array(load("num_list2.mat","num_list2"));
% 
% 
% illumination_na = sqrt(kx_relative.^2+ky_relative.^2);
% na_used2=illumination_na(num_list2);
% %column_sort=sort(na_used);
% na_used1 = na_used2(:,1)+na_used2(:,2);
% %head_one=column_sort(1,:);
% [dis_lit_new,idx_led_new] = sort(na_used1);
% % Nsh_lit = zeros(numlit,Nimg);
% % Nsv_lit = zeros(numlit,Nimg);
% 
% 
% num_list2 = num_list2(idx_led_new,:);

num_list2 = zeros(num_image,2);
look_1_9 = zeros(num_image,2);
bright = [127,128,129,146,147,148,165,166,167];
newimSeqLowRes2 = zeros(128,128,num_image);
k=1;
while num_list2(num_image,1)==0
    i = randi([1,9],1,1,'int8');
    j = randi([1,9],1,1,'int8');

    if i~=j 
            num_i = bright(i);
            num_j = bright(j);
            num_list2(k,1) = num_i;
            num_list2(k,2) = num_j;
            look_1_9(k,1) = i;
            look_1_9(k,2) = j;
            newimSeqLowRes2(:,:,k) = newimSeqLowRes(:,:,i)+newimSeqLowRes(:,:,j);
            k = k+1;
    end
end 

%objectlow = newimSeqLowRes(:,:,1);
%I_one =(newimSeqLowRes(:,:,147));
objectlow = sqrt(abs(newimSeqLowRes2(:,:,1)));
figure(85);
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
%seq293 = struct2array(load('seq288.mat', 'seq293'));%����������
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


    %% rest led zheng8
for tt=1:loop
    
    iter = iter+1;
    err2 = 0; 
    
    
    for i3 = 1:num_image
%         i2=seq293(i3);
% 
%         [x,y]=find(num_list==i2);
%         if mod(i3,8) ==0
%            gamma2 = gamma2/2;
%         end


%%%%%
        lowResFT2 = zeros(m1,m1,2);
        im_lowRes2 = zeros(m1,m1,2);
        maxP = zeros(n,m);
        for r=1:2
            kxc=round((n+1)/2+kx(1,num_list2(i3,r))/dkx);%round�Ƕ������ڵ�Ԫ��ȡ��
            kyc=round((m+1)/2+ky(1,num_list2(i3,r))/dky);
            kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
            kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
            kxc2=round((n+1)/2-kx(1,num_list2(i3,r))/dkx);%round�Ƕ������ڵ�Ԫ��ȡ��
            kyc2=round((m+1)/2-ky(1,num_list2(i3,r))/dky);
            kyl2=round(kyc2-(m1-1)/2);kyh2=round(kyc2+(m1-1)/2);
            kxl2=round(kxc2-(n1-1)/2);kxh2=round(kxc2+(n1-1)/2);
            maxP(kyl:kyh,kxl:kxh) = pupil;
            pupil2=maxP(128:128+128-1,128:128+128-1);
            lowResFT2_2(:,:,r)=objectRecoverFT(128:128+128-1,128:128+128-1).*pupil2.*pupil0;
            lowResFT2(:,:,r)=objectRecoverFT(kyl2:kyh2,kxl2:kxh2).*pupil.*pupil0;
            im_lowRes2(:,:,r) = (ifft2(ifftshift(lowResFT2(:,:,r))));
            im_lowRes2_2(:,:,r) = (ifft2(ifftshift(lowResFT2_2(:,:,r))));
        end
                
        %im_lowRes_pinpu=lowResFT_1+lowResFT_2+lowResFT_3;
        im_lowRes_oI=(abs(im_lowRes2(:,:,1))).^2+(abs(im_lowRes2(:,:,2))).^2;
        im_lowRes_oI_2=(abs(im_lowRes2_2(:,:,1))).^2+(abs(im_lowRes2_2(:,:,2))).^2;
        
%         kxc=round((n+1)/2+kx(1,i2)/dkx);%round�Ƕ������ڵ�Ԫ��ȡ��
%         kyc=round((m+1)/2+ky(1,i2)/dky);
%         kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
%         kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
%         lowResFT_ori=(m1/m)^2*objectRecoverFT(kyl:kyh,kxl:kxh).*pupil;
%         im_lowRes_ori = ifft2(ifftshift(lowResFT_ori));
                       
        %im_lowRes = (m/m1)^2*imSeqLowRes(:,:,i2).*exp(1i.*angle(im_lowRes));%�滻���
        %im_lowRes = (m/m1)^2*sqrt(abs(newimSeqLowRes(:,:,x))./max(max(abs(newimSeqLowRes(:,:,x)))))./sqrt(abs(im_lowRes_oI)./max(max(abs(im_lowRes_oI)))).*im_lowRes_ori;%�滻���
        im_lowResk = zeros(m1,m1,2);
        for kk=1:2
             im_lowResk(:,:,kk) = sqrt(abs((newimSeqLowRes2(:,:,i3))))./sqrt(abs(im_lowRes_oI)).*im_lowRes2(:,:,kk).*pupil.*pupil0;%�滻���
             im_lowResk_2(:,:,kk) = sqrt(abs((newimSeqLowRes2(:,:,i3))))./sqrt(abs(im_lowRes_oI_2)).*im_lowRes2_2(:,:,kk).*pupil.*pupil0;
        end
        lowResFT=fftshift(fft2(im_lowResk));%ת��Ϊ����Ҷ��
        lowResFT_2=fftshift(fft2(im_lowResk_2));
        OP_diff=(m1/m)^2*lowResFT-lowResFT2;
        OP_diff2=(m1/m)^2*lowResFT_2-lowResFT2_2;
        dobject = zeros(n,m);
        dotest = zeros(n,m);
        dpupil = 0;
        sumobject = 0;
        sumpupil = zeros(n,m);
        for kk=1:2
            kxc=round((n+1)/2+kx(1,num_list2(i3,kk))/dkx);%round�Ƕ������ڵ�Ԫ��ȡ��
            kyc=round((m+1)/2+ky(1,num_list2(i3,kk))/dky);
            kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
            kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
            kxc2=round((n+1)/2-kx(1,num_list2(i3,r))/dkx);%round�Ƕ������ڵ�Ԫ��ȡ��
            kyc2=round((m+1)/2-ky(1,num_list2(i3,r))/dky);
            kyl2=round(kyc2-(m1-1)/2);kyh2=round(kyc2+(m1-1)/2);
            kxl2=round(kxc2-(n1-1)/2);kxh2=round(kxc2+(n1-1)/2);
            
            OFTdown = objectRecoverFT(kyl2:kyh2,kxl2:kxh2);
            sumobject = sumobject+abs(OFTdown).^2;
            
            maxP = zeros(n,m);
      
            maxP(kyl:kyh,kxl:kxh) = maxP(kyl:kyh,kxl:kxh)+pupil;
            pupil2=maxP(128:128+128-1,128:128+128-1);
            sumpupil(kyl:kyh,kxl:kxh) = sumpupil(kyl:kyh,kxl:kxh)+abs(pupil).^2;
            
            
            %dobject=dobject+abs(pupil2).*conj(pupil2).*(OP_diff2(:,:,kk));
            dobject(128:128+128-1,128:128+128-1)=dobject(128:128+128-1,128:128+128-1)+abs(pupil2).*conj(pupil2).*(OP_diff2(:,:,kk));
     
            dpupil =  dpupil+(abs(OFTdown).*conj(OFTdown)).*(OP_diff(:,:,kk));

        end
        
        pupil =(pupil +gamma3*dpupil/max(max(abs(objectRecoverFT)))./(sumobject+beta)).*pupil0;%(abs(F((abs(Ft(Ps)).^2)))./max(max(abs(F((abs(Ft(Ps)).^2))))));%1,9e14;
        
        %pupil =(pupil+dptest).*Ps;
        %dotest=.1*(dO./(max(max(abs(P)))))./(sumP+alpha);%0.8,80/384/384;
        objectRecoverFT = objectRecoverFT + gamma2*(dobject./max(max(abs(pupil))))./(sumpupil+alpha);
        
        
        %pupil =  pupil+gamma3/max(max(abs(objectRecoverFT))).*(abs(OFTdown).*conj(OFTdown)).*(OP_diff(:,:,kk))./(abs(OFTdown).^2+beta).*pupil0;

           
        err2 = err2+sqrt(sum(sum(((abs(newimSeqLowRes2(:,:,i3))-abs(im_lowRes_oI)).^2))));

    end

%

    fprintf('| %2d   | %.2e |\n',iter,err2);
    objectRecover=ifft2(ifftshift(objectRecoverFT));
    figure(51);
    subplot(221); imagesc(abs(objectRecover)); axis image; colormap gray; colorbar;
    title('ampl(o)');
    subplot(222); imagesc(angle(objectRecover)); axis image; colormap gray; colorbar;
    title('phase(o)');
    subplot(223); imagesc(abs(log((objectRecoverFT+1)))); axis image; colormap gray; colorbar;
    title('objectRecoverFT');
    subplot(224); imagesc(abs(log(objectRecoverFT0+1))); axis image; colormap gray; colorbar;
    title('objectRecoverFT0');%pupil    log(objectRecoverFT0)

end
objectRecover=ifft2(ifftshift(objectRecoverFT));
figure(52);
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
%imwrite(ampl,'/home/jnu733/qwm/result/matlab_result/2led/usaf2ampli_36_tian.png','png');
% fileName = ['simulate293ampli.tif'];
% fullPath = fullfile(savePath,fileName); % ��ȡ�����ı���·�����ļ���
% imwrite(I_est,fullPath, 'tif', 'Compression', 'none'); % ����Ϊ16λ��ͨ��tifͼƬ���ر�ѹ��

ampl=angle(objectRecover);
am_high = max(max(ampl));
am_low = min(min(ampl));
ampl = (ampl-am_low)./(am_high-am_low);
%imwrite(ampl,'/home/jnu733/qwm/result/matlab_result/2led/usaf2angle_36_tian.png','png');
% fileName = ['simulate293angle.tif'];
% fullPath = fullfile(savePath,fileName); % ��ȡ�����ı���·�����ļ���
% imwrite(I_est,fullPath, 'tif', 'Compression', 'none'); % ����Ϊ16λ��ͨ��tifͼƬ���ر�ѹ��



