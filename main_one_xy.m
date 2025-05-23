%�����һ���Ƶĸ���Ҷ���㣬��kx,ky����Ns
maxIter = 30;
alpha = 0.01;%1
beta = 1e20;%1e3

%% Reconstruction library locates here
%clear all;
addpath(['./FP_Func']);
addpath('natsortfiles');

filedir = ['/home/jnu733/qwm/data/usaf1led/'];
% filedir = ['D:/workapp/matlab/workstation/QWM-FPM/data/1_LED/tif/'];
imglist = dir([filedir,'*.tif']);
%nstart = [1017,1217];
nstart = [973,1173];
%nstart = [800,1400];
N = natsortfiles({imglist.name});
numlit = 1;
n1 = 2160; n2 = 2560;
%% read in all images into the memory first
fprintf(['loading the images...\n']);
tic;
Nimg = length(imglist);
Iall = zeros(n1,n2,Nimg);
Ibk = zeros(Nimg,1);
for m = 1:Nimg
    fn = [filedir,N{m}];  %added for sorting purpose
    disp(fn);
    % all image data
    I = double(imread(fn));
    Iall(:,:,m) = I  ;
    bk1 = mean2(double(Iall(140:200,140:200,m)));
    bk2 = mean2(double(Iall(150:210,2420:2480,m)));%(492:600,2380:2520,m)
    Ibk(m) = mean([bk1,bk2]);
    if Ibk(m)>1000 
        Ibk(m) = Ibk(m-1);
    end
    
end
    
fprintf(['\nfinish loading images\n']);
toc;

%% define processing ROI
Np = 128;

F = @(x) fftshift(fft2(x));
Ft = @(x) ifft2(ifftshift(x));
row = @(x) x(:).';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelength of illumination, assume monochromatic
% R: 624.4nm +- 50nm
% G: 518.0nm +- 50nm
% B: 476.4nm +- 50nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 0.6292;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical aperture of the objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NA = 0.1;
% maximum spatial frequency set by NA
um_m = NA/lambda;
% system resolution based on the NA
dx0 = 1/um_m/2;
mag = 8.1485;

dpix_c = 6.5; %6.5um pixel size on the sensor plane
dpix_m = dpix_c/mag; 
% FoV in the object space
FoV = Np*dpix_m;
% sampling size at Fourier plane set by the image size (FoV)
% sampling size at Fourier plane is always = 1/FoV
if mod(Np,2) == 1
    du = 1/dpix_m/(Np-1);
else
    du = 1/FoV;
end

m = 1:Np;
[mm,nn] = meshgrid(m-round((Np+1)/2));
ridx = sqrt(mm.^2+nn.^2);
um_idx = um_m/du;
% assume a circular pupil function, lpf due to finite NA
w_NA = double(ridx<um_idx);
pupil = w_NA;

clear m mm nn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up image corrdinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncent = [1080,1280];
img_ncent = nstart-ncent+Np/2;
img_center = (nstart-ncent+Np/2)*dpix_m;
img_start = nstart*dpix_m;
img_end = (nstart+Np)*dpix_m;

%% LED array geometries and derived quantities

ds_led = 4e3; %4mm
z_led = 67.5e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diameter of # of LEDs used in the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dia_led = 19;
lit_cenv = 13;
lit_cenh = 14;
vled = [0:31]-lit_cenv;
hled = [0:31]-lit_cenh;

[hhled,vvled] = meshgrid(hled,vled);
rrled = sqrt(hhled.^2+vvled.^2);
LitCoord = rrled<dia_led/2;

% total number of LEDs used in the experiment
Nled = sum(LitCoord(:));
% index of LEDs used in the experiment
Litidx = find(LitCoord);
[xmid,ymid] = find(LitCoord);
dd2 = sqrt((-(xmid-14)*ds_led-img_center(1)).^2+(-(ymid-15)*ds_led-img_center(2)).^2+z_led.^2);
kx = round((-(xmid-14)*ds_led-img_center(1))./dd2/lambda/du);
ky = round((-(ymid-15)*ds_led-img_center(2))./dd2/lambda/du);

% corresponding angles for each LEDs
dd = sqrt((-hhled*ds_led-img_center(1)).^2+(-vvled*ds_led-img_center(2)).^2+z_led.^2);
sin_thetav = (-hhled*ds_led-img_center(1))./dd;
sin_thetah = (-vvled*ds_led-img_center(2))./dd;

illumination_na = sqrt(sin_thetav.^2+sin_thetah.^2);
% corresponding spatial freq for each LEDs
%
vled = sin_thetav/lambda;
uled = sin_thetah/lambda;
% spatial freq index for each plane wave relative to the center
idx_u = round(uled/du);
idx_v = round(vled/du);

illumination_na_used = illumination_na(LitCoord);

% number of brightfield image
NBF = sum(illumination_na_used<NA);
% maxium spatial frequency achievable based on the maximum illumination
% angle from the LED array and NA of the objective
um_p = max(illumination_na_used)/lambda+um_m;
% resolution achieved after freq post-processing
dx0_p = 1/um_p/2;

disp(['synthetic NA is ',num2str(um_p*lambda)]);
% assume the max spatial freq of the original object
% um_obj>um_p
% assume the # of pixels of the original object
N_obj = round(2*um_p/du)*2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to enforce N_obj/Np = integer to ensure no FT artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_obj = ceil(N_obj/Np)*Np;

Imea = double(Iall(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1,:));

ledidx = 1:Nled;
lit = Litidx(ledidx);%
lit = reshape(lit,1,293);
[dis_lit2,idx_led] = sort(reshape(illumination_na_used,1,Nled));%sort by NA

Nsh_lit = zeros(numlit,Nimg);
Nsv_lit = zeros(numlit,Nimg);
for m = 1:Nimg
    % corresponding index of spatial freq for the LEDs are lit
    lit0 = lit(1,m);
    Nsh_lit(1,m) = idx_u(lit0);
    Nsv_lit(1,m) = idx_v(lit0);
end

% reorder the LED indices and intensity measurements according the previous
% dis_lit
Ns = [];
Ns(1,:,1) = Nsv_lit;
Ns(1,:,2) = Nsh_lit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing the data to DENOISING is IMPORTANT
% background subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = Imea(:,:,idx_led);
Ibk = Ibk(idx_led);
Ns = Ns(:,idx_led,:);
kx = kx(idx_led);
ky = ky(idx_led);

for m = 1:Nimg
    Itmp = I(:,:,m);
    Itmp = Itmp-Ibk(m);
%     Itmp = awgn(Itmp,0,'measured');
    Itmp(Itmp<0) = 0;
    I(:,:,m) = Itmp;
end

I_show = I(:,:,1);

tol = 1;%1     
minIter = 2;
monotone = 1;  %
display = 'full';%0;%'iter';
upsamp = @(x) padarray(x,[(N_obj-Np)/2,(N_obj-Np)/2]);
O = F(sqrt(I(:,:,1)));
O = upsamp(O);
Fbegin =O;
P = w_NA;
Ps = w_NA;
StepSize = 0.1;

T0 = clock;

fprintf('| iter |  rmse    |\n');
for j=1:20, fprintf('-'); end
fprintf('\n');

err1 = inf;
err2 = 50;
err = [];
iter = 0;

if strcmp(display,'full')
    
    o = Ft(O);
    f1 = figure(87);
    subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
    title('ampl(o)');
    subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
    title('phase(o)');
    subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
    title('ampl(P)');
    subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
    title('phase(P)');
    drawnow;
end

%% main algorithm starts here
fprintf('| %2d   | %.2e |\n',iter,err1);

% operator to crop region of O from proper location at the O plane
downsamp = @(x,cen) x(cen(1)-floor(Np/2):cen(1)-floor(Np/2)+Np-1,...
    cen(2)-floor(Np/2):cen(2)-floor(Np/2)+Np-1);

while iter<maxIter
    err1 = err2;
    err2 = 0;
    iter = iter+1;
    cen0 = round((N_obj+1)/2);
    for m = 1:Nimg
        % initilize psi for correponing image, ROI determined by cen
        Psi0 = zeros(Np,Np); %128x128
        cen = zeros(2,1);        
        %cen(:,1) = cen0-row(Ns(1,m,:));
        cen(1,1) = cen0-ky(m);
        cen(2,1) = cen0-kx(m);
        Psi0(:,:) = downsamp(O,cen(:,1)).*P;%FF P=128*128
        
        I_mea = I(:,:,m);
        psi0 = Ft(Psi0);%%RR 128x128
        I_est = abs(psi0).^2; %caculate the error
        Psi = zeros(Np,Np);%128*128
        
        Psi = F(sqrt(I_mea).*exp(1i .* angle(psi0)));
        %psi(:,:) = F(sqrt(I_mea).*psi0(:,:)./sqrt(I_est+eps)); %%
        
        dpsi = Psi-Psi0;%128*128
        %Omax = abs(O(cen0(1),cen0(2)));
        Omax = max(max(O));

        % operator to put P at proper location at the O plane
        n1 = cen-floor(Np/2);
        n2 = n1+Np-1;
        % operator to crop region of O from proper location at the O plane
        %downsamp2 = @(x) x(n1(1):n2(1),n1(2):n2(2));

        O1 = O(n1(1):n2(1),n1(2):n2(2));

        O(n1(1):n2(1),n1(2):n2(2)) = O(n1(1):n2(1),n1(2):n2(2))...
            + 0.1 * 1/max(max(abs(P)))*abs(P).*conj(P).*dpsi./(abs(P).^2+alpha);
        P = P+0.5/Omax*(abs(O1).*conj(O1)).*dpsi./(abs(O1).^2+beta).*Ps;
        err2 = err2+sqrt(sum(sum((I_mea-I_est).^2)));
    end
    err = [err,err2];

    if strcmp(display,'full')
        o = Ft(O);
        figure(89);
        subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
        title('ampl(o)');
        subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
        title('phase(o)');
        subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
        title('ampl(P)');
        subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
        title('phase(P)');
        drawnow;
    end
    Fend = O;
    figure(99);
    subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
    title('ampl(o)');
    subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
    title('phase(o)');
    subplot(223); imagesc(abs(log(Fbegin+1))); axis image; colormap gray; colorbar;
    title('Fbegin');
    subplot(224); imagesc(abs(log(Fend+1))); axis image; colorbar;
    title('Fend');
    drawnow;
    
    
    fprintf('| %2d   | %.2e |\n',iter,err2);
end
fprintf('processing completes\n');
 
ampl = abs(o);
high = max(max(ampl));
low = min(min(ampl));
ampl = (ampl-low)./(high-low);
imwrite(ampl,'result.png');
