%% tianlei 1led liang && 4led anchang
%clc;
clc;
clear all;
%close all;
addpath(['./FP_Func']);
addpath('natsortfiles');

loop=50;
alpha = eps;%1
beta = eps;%1e16;%1e20;%1e3
%gamma1 =8;
gamma2 =.0001 ;%object usaf 0.01
gamma3=.001;%pupil usaf 0.1

LED_num=2;%8;
Img_assemble_num=290;

LEDheight=67.5e3;
LEDgap=4e3;
waveLength = 0.6292;
m = 384; 
n = 384;

m1 = m/3; %
n1 = n/3; %

mag = 8.1485;
dpix_c = 6.5; %6.5um pixel size on the sensor plane
dpix_m = dpix_c/mag; 
dkx = 1/((m1-1)*dpix_m);
dky = 1/((m1-1)*dpix_m);
NA = 0.1;
% maximum spatial frequency set by NA
um_m = NA/waveLength;
rx = 1:m1;
[rxx,ryy] = meshgrid(rx-round((m1+1)/2));
ridx = sqrt(rxx.^2+ryy.^2);

rx2 = 1:m;
[rxx2,ryy2] = meshgrid(rx2-round((m+1)/2));
ridx2 = sqrt(rxx2.^2+ryy2.^2);
um_idx = um_m/dkx;
% assume a circular pupil function, lpf due to finite NA
CTF = double(ridx<um_idx);
CTF2 = double(ridx2<um_idx);

%nstart = [973,1173];
%nstart = [888,1088];
nstart=[800,1400];
%nstart = [1017,1217];
ncent = [1080,1280];
%% read in all images into the memory first
%filedir1 = 'F:/tif2/';
%filedir1='./8LED/tif/';
filedir1 = '/home/jnu733/qwm/data/cell1led/';
%filedir2 = '/home/jnu733/qwm/data/cell1led/';

% D:/workapp/matlab/workstation/QWM-FPM/data/8led/tif/
% /home/jnu733/qwm/data
imglist1 = dir([filedir1,'*.tif']);
N = natsortfiles({imglist1.name});
Nimg = length(imglist1);
%img_ncent = nstart-ncent+m1/2;
img_center = (nstart-ncent+m1/2)*dpix_m;
xlo = struct2array(load('xloylo.mat', 'xlo'));%
ylo = struct2array(load('xloylo.mat', 'ylo'));%
load('./expt_lit-8.mat');




%numlit = 8;
%nn1 = 128; nn2 = 128;

fprintf('loading the images...\n');
tic;

Iall = zeros(2160,2560,Nimg,'uint16');
Iall2 = zeros(2160,2560,Nimg,'uint16');

%[a,mm3]=sort(k_all);
mm=[147,148,129,128,127,146,165,166,167,168,149,130,111,110,109,108,107,126,145,164,183,184,185,186,187,188, ...
    169,150,131,112,93,92,91,90,89,88,87,106,125,144,163,182,201,202,203,204,205,206,207,208,189,170,151,132,113,94, ...
    76,75,74,73,72,71,70,69,68,86,105,124,143,162,181,200,218,219,220,221,222,223,224,225,226,227,209,190,171,152,133,114,95,77,60, ...
    59,58,57,56,55,54,53,52,51,50,67,85,104,123,142,161,180,199,217,234,235,236,237,238,239,240,241,242,243,244,245, ...
    228,210,191,172,153,134,115,96,78,61,45,44,43,42,41,40,39,38,37,36,35,34,33,49,66,84,103,122,141,160,179,198,216,233,249, ...
    250,251,252,253,254,255,256,257,258,259,260,261,262,246,229,211,192,173,154,135,116,97,79,62,46,31,30,29,28,27,26,25,24,23,22,21,20,19, ...
    32,48,65,83,102,121,140,159,178,197,215,232,248,263,264,265,266,267,268,269,270,271,272,273,274,275,247,230,212,193,174,155,136,117,98,80,63, ...
    18,17,16,15,14,13,12,11,10,9,8,47,64,82,101,120,139,158,177,196,214,231, ...
    276,277,278,279,280,281,282,283,284,285,286,213,194,175,156,137,118,99,7,6,5,4,3,2,1, ...
    81,100,119,138,157,176,195,287,288,289,290,291,292,293];



Ibk = zeros(Nimg,1);
for mm1 = 1:Nimg
    fn = [filedir1,N{mm1}];  %added for sorting purpose
    disp(fn);
    I = double(imread(fn));

    bk1 = mean2(double(I(1:100,1:100)));
    bk2= mean2(double(I(492:600,2380:2520)));
    Ibk(mm1) = mean([bk1,bk2]);
    
    
    Iall2(:,:,mm1) = I  ;

       
      if Ibk(mm1)>300 && mm1 >2
             Ibk(mm1) = Ibk(mm1-1);
      end
%      if Ibk(mm1)>1000&& mm1 >2
%             Ibk(mm1) = Ibk(mm1-1);
%      end
%      
     
end

%Iall2=zeros(2160,2560,Nimg,'uint16');
%Iall2(:,:,mm(1))=Iall(:,:,mm(1));
%for r = 2:Nimg
%    Iall2(:,:,mm(r))=mean2(Iall(:,:,mm(r-1)))/mean2(Iall(:,:,mm(r))).*Iall(:,:,mm(r));
%end 

% Ibk(1) = Ibk(3);
% Ibk(2) = Ibk(3);
%Ibk(:)=800;
for mm1 = 1:Nimg
    Itmp = Iall2(:,:,mm1);
    Itmp = Itmp-Ibk(mm1);
%     Itmp = awgn(Itmp,0,'measured');
    Itmp(Itmp<0) = 0;
    Iall2(:,:,mm1) = Itmp;
    
end

Iallend = double(Iall2(nstart(1):nstart(1)+m1-1,nstart(2):nstart(2)+m1-1,:));
%clear I Iall


I_resize = double(Iallend);


objectlow_first = sqrt(abs(I_resize(:,:,146)));
figure(85);
subplot(221); imagesc(abs(objectlow_first)); axis image; colormap gray; colorbar;
title('ampl(o)');
subplot(222); imagesc(angle(objectlow_first)); axis image; colormap gray; colorbar;
title('phase(o)');
subplot(223); imagesc(angle(objectlow_first)); axis image; colormap gray; colorbar;
title('ampl(O)');
subplot(224); imagesc(angle(objectlow_first)); axis image; colorbar;
title('phase(O)');


F = @(x) fftshift(fft2(x));
Ft = @(x) ifft2(ifftshift(x));
upsamp = @(x) padarray(x,[(m-m1)/2,(m-m1)/2]);




%LED8293=reshape(ledidx,8,293);
%LEDx=reshape(xlit,8,293);
%LEDy=reshape(ylit,8,293);
%dd =sqrt(((LEDx-13).*LEDgap-img_center(1)).^2+((LEDy-14).*LEDgap-img_center(1)).^2+LEDheight^2);
%kxx=((LEDx-13)*LEDgap-img_center(1))./dd/waveLength;
%kyy=((LEDy-14)*LEDgap-img_center(1))./dd/waveLength;
%kall=kxx.^2+kyy.^2;
%kall_sort_1=sort(kall,1);
%[LED_list1,LED_list2]=sort(kall_sort_1(1,:),2);
%[LED_list1,LED_list2]=sort(sum(kall_sort_1,1));

LEDij=zeros(LED_num,Img_assemble_num);
I_assemble=zeros(n1,m1,Img_assemble_num);
dd_assemble=zeros(LED_num,Img_assemble_num);
kxx=zeros(LED_num,Img_assemble_num);
kyy=zeros(LED_num,Img_assemble_num);

for ii=1:Img_assemble_num
    jj=ii+10;%;randi([1,15],1,1);
    if jj>293
        jj=jj-14;
    end
    I_assemble(:,:,ii) = I_resize(:,:,mm(ii))+I_resize(:,:,mm(jj));
    LEDij(1,ii)=mm(ii);
    LEDij(2,ii)=mm(jj);
    i1=find(ledidx==mm(ii)); 
    i2=find(ledidx==mm(jj));
    dd_assemble(1,ii) =sqrt(((xlit(i1(1))-13)*LEDgap-img_center(1))^2+((ylit(i1(1))-14)*LEDgap-img_center(2))^2+LEDheight^2);
    dd_assemble(2,ii) =sqrt(((xlit(i2(1))-13)*LEDgap-img_center(1))^2+((ylit(i2(1))-14)*LEDgap-img_center(2))^2+LEDheight^2);
    kxx(1,ii)=((ylit(i1(1))-14)*LEDgap-img_center(2))./dd_assemble(1,ii)/waveLength;
    kyy(1,ii)=((xlit(i1(1))-13)*LEDgap-img_center(1))./dd_assemble(1,ii)/waveLength;
    kxx(2,ii)=((ylit(i2(1))-14)*LEDgap-img_center(2))./dd_assemble(2,ii)/waveLength;
    kyy(2,ii)=((xlit(i2(1))-13)*LEDgap-img_center(1))./dd_assemble(2,ii)/waveLength;

end 

%load('./img_re.mat');


I_RecoverFT2 = Ft(upsamp(F(sqrt(I_assemble(:,:,1)))));
phase=abs(I_RecoverFT2);

ymax=0.5;ymin=-0.5;
xmax = max(max(phase)); %求得InImg中的最大值
xmin = min(min(phase)); %求得InImg中的最小值
phase = (ymax-ymin)*(phase-xmin)/(xmax-xmin) + ymin;
I_RecoverFT = abs(I_RecoverFT2).*exp(1i.* phase);
I_RecoverFT =F(I_RecoverFT );

pupil0=double(CTF);
pupil=pupil0;
maxP = CTF2;

fprintf('| iter |  rmse    |\n');
for j=1:20, fprintf('-'); end
fprintf('\n');
err2 = 0;
err = [];
iter = 0;

fprintf('| %2d   | %.2e |\n',iter,err2);

for tt=1:loop
    iter = iter+1;
    err2 = 0; 

    for ii = 1:Img_assemble_num
        kxc=round((n+1)/2+kxx(:,ii)./dkx);
        kyc=round((m+1)/2+kyy(:,ii)./dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
 
        lowResFT = zeros(m1,m1,LED_num);
        img_lowRes = zeros(m1,m1,LED_num);
        img_lowRes_sum=zeros(m1,m1);
       
        for r=1:LED_num
            
            lowResFT(:,:,r)=I_RecoverFT(kyl(r,:):kyh(r,:),kxl(r,:):kxh(r,:)).*pupil;
            img_lowRes(:,:,r) = (ifft2(ifftshift(lowResFT(:,:,r))));
            img_lowRes_sum=img_lowRes_sum+(abs(img_lowRes(:,:,r))).^2;
        end
        img_lowRes_reset = zeros(m1,m1,LED_num);
        
        for r=1:LED_num
        
        img_lowRes_reset(:,:,r) = sqrt(abs(I_assemble(:,:,ii)))./sqrt(abs(img_lowRes_sum(:,:))).*img_lowRes(:,:,r);
        end
        lowRes_resetFT=fftshift(fft2(img_lowRes_reset));

        OP_diff=lowRes_resetFT-lowResFT;
        
        dobject = zeros(n,m);
        dotest = zeros(n,m);
        dpupil = zeros(m1,m1);
        sumobject = 0;
        sumpupil = zeros(n,m);
        for kk=1:LED_num
            
            
            OFTdown = I_RecoverFT(kyl(kk,:):kyh(kk,:),kxl(kk,:):kxh(kk,:));
            sumobject = sumobject+abs(OFTdown).^2;

            sumpupil(kyl(kk,:):kyh(kk,:),kxl(kk,:):kxh(kk,:)) = sumpupil(kyl(kk,:):kyh(kk,:),kxl(kk,:):kxh(kk,:))+abs(pupil).^2;
           
            dobject(kyl(kk,:):kyh(kk,:),kxl(kk,:):kxh(kk,:))=dobject(kyl(kk,:):kyh(kk,:),kxl(kk,:):kxh(kk,:))+abs(pupil).*conj(pupil).*(OP_diff(:,:,kk));
     
            dpupil =  dpupil+(abs(OFTdown).*conj(OFTdown)).*(OP_diff(:,:,kk));

        end
        
        pupil =(pupil +gamma3*dpupil/max(max(abs(I_RecoverFT)))./(sumobject+beta)).*pupil0;
        I_RecoverFT = I_RecoverFT + gamma2*(dobject./max(max(abs(pupil))))./(sumpupil+alpha);
        
        err2 = err2+sqrt(sum(sum(((abs(I_assemble(:,:,ii))-abs(img_lowRes_sum)).^2))));

    end

%

    fprintf('| %2d   | %.2e |\n',iter,err2);
    objectRecover=ifft2(ifftshift(I_RecoverFT));
    figure(51);
    subplot(221); imagesc(abs(objectRecover)); axis image; colormap gray; colorbar;
    title('ampl(o)');
    subplot(222); imagesc(angle(objectRecover)); axis image; colormap gray; colorbar;
    title('phase(o)');
    subplot(223); imagesc(log(abs(dobject(:,:,1))+eps)); axis image; colormap gray; colorbar;
    title('objectRecoverFT');
    subplot(224);  imagesc((angle(OP_diff(:,:,2))+eps)); axis image; colormap gray; colorbar;
    title('objectRecoverFT0');%pupil    log(objectRecoverFT0)


end
objectRecover=ifft2(ifftshift(I_RecoverFT));
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

