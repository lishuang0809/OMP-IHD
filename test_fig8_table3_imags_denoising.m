clc; clear; close all;

rng(30);  % 20

%% Load images
img0 = imread('images/MRI3_k16_True.jpg');  % Read the image
% Check if the image is RGB (3D) and convert to grayscale if necessary
if size(img0, 3) == 3
    img0 = rgb2gray(img0);
end

img_clean = im2double(img0); % Normalize to [0,1]
[N0,J0] = size(img_clean);

% Patch Parameters
patch_size = 8;
n_patches = (N0 / patch_size) * (J0 / patch_size); % # of patches

N = patch_size^2; J = n_patches;
k = 16;

% Generate Overcomplete DCT Dictionary (64x676)
RR = 2; % redundancy factor used for setting dictionary's size
DCT_dict = TwoD_DCT_dic_BaiH(patch_size, RR);
L = size(DCT_dict,2);
Q = DCT_dict;

X_clean = zeros(N,J);
count = 1;
for i = 1:patch_size:N0
    for j = 1:patch_size:J0
        patch = img_clean(i:i+patch_size-1, j:j+patch_size-1); % Extract 8x8 patch
        X_clean(:, count) = patch(:); % Flatten and store as column
        count = count + 1;
    end
end

% Add low-dimensional noise to image
sigma_e = 0.2;
N_e = 20;
Gamma = randn(N);
[Gamma,~,~]=svd(Gamma);
Gamma = Gamma(:,1:N_e);
V_alpha = randn(N_e,J); 
V_alpha = normc(V_alpha);
V_alpha = V_alpha*sigma_e;
E = Gamma*V_alpha;          % Note ||E(:,j)||_2 = sigma_e

X = X_clean + E;

% Reconstruct the image by placing patches back
img_noise = zeros(size(img0));
count = 1;
for i = 1:patch_size:N0
    for j = 1:patch_size:J0
        patch = reshape(X(:, count), [patch_size, patch_size]); % Convert column to 8x8 patch
        img_noise(i:i+patch_size-1, j:j+patch_size-1) = patch; % Place back in correct position
        count = count + 1;
    end
end

I_noisy = im2uint8(img_noise);

psnr_noi = psnr(I_noisy, img0);

figure(1); clf;
fs = 20;
subplot(1,2,1)
imshow(img0);
title('Original','FontSize',fs)

subplot(1,2,2)
imshow(I_noisy);
title(['Noisy (psnr = ', num2str(psnr_noi) ')'], 'FontSize', fs);

%% Preparation 

[U,Sig,V]=svd(Q);  Sig=Sig(:,1:N);        % Q = U[Sig 0]V^T
V_1 = V(:,1:N);                % Q_bar = V_1*V_1';
T_0 = V(:,1:N)/Sig*U';         %SL: was V(:,1:N)*inv(Sig)*U';
W = V(:,N+1:L);
W_c = Sig\U';                  %SL: was inv(Sig)*U';
S_0 = T_0*X;

if sigma_e==0
    W_bar = W;
else
    W_bar = [W -T_0*Gamma];
end

%% Sparse recovery with different algorithms

% ======= Yin_IRLS ========
alg = 104;
[S_Yin_IRLS,time_Yin_IRLS] = NewAlgs_4_Noise_imag(Q,X,k,W_bar,S_0,V_1,alg);
disp(['Time_IRLS: ', num2str(time_Yin_IRLS), ' seconds']);

% ======= Alg_{GLQ} ========
alg = 109;
[S_GLQ,time_GLQ] = NewAlgs_4_Noise_imag(Q,X,k,W_bar,S_0,V_1,alg);
disp(['Time_GLQ: ', num2str(time_GLQ), ' seconds']);


%%
X_IRLS0 = Q*S_Yin_IRLS;
X_GLQ0 = Q*S_GLQ;


% Reconstruct the image by placing patches back
X_IRLS = zeros(size(img0));
X_GLQ = zeros(size(img0));
count = 1;
for i = 1:patch_size:N0
    for j = 1:patch_size:J0
        % patch = reshape(X_OMP0(:, count), [patch_size, patch_size]); % Convert column to 8x8 patch
        % X_OMP(i:i+patch_size-1, j:j+patch_size-1) = patch; % Place back in correct position
        patch = reshape(X_IRLS0(:, count), [patch_size, patch_size]); % Convert column to 8x8 patch
        X_IRLS(i:i+patch_size-1, j:j+patch_size-1) = patch; % Place back in correct position
        % patch = reshape(X_GL20(:, count), [patch_size, patch_size]); % Convert column to 8x8 patch
        % X_GL2(i:i+patch_size-1, j:j+patch_size-1) = patch; % Place back in correct position
        patch = reshape(X_GLQ0(:, count), [patch_size, patch_size]); % Convert column to 8x8 patch
        X_GLQ(i:i+patch_size-1, j:j+patch_size-1) = patch; % Place back in correct position
        % patch = reshape(X_FGLQ0(:, count), [patch_size, patch_size]); % Convert column to 8x8 patch
        % X_FGLQ(i:i+patch_size-1, j:j+patch_size-1) = patch; % Place back in correct position        
        count = count + 1;
    end
end

I_IRLS = im2uint8(X_IRLS);
I_GLQ = im2uint8(X_GLQ);

psnr_IRLS = psnr(I_IRLS, img0);
psnr_GLQ = psnr(I_GLQ, img0);


%%
figure(2); clf;
fs = 20;
subplot(2,4,1)
imshow(img0);
title('Original','FontSize',fs)

subplot(2,4,2)
imshow(I_noisy);
title(['Noisy (psnr = ', num2str(psnr_noi) ')'], 'FontSize', fs);

subplot(2,4,3)
imshow(I_OMP);
title(['OMP (psnr = ', num2str(psnr_OMP) ')'], 'FontSize', fs);
xlabel(['Time = ', num2str(time_OMP)], 'FontSize', fs);

subplot(2,4,4)
imshow(I_IRLS);
title(['IRLS (psnr = ', num2str(psnr_IRLS) ')'], 'FontSize', fs);
xlabel(['Time = ', num2str(time_Yin_IRLS)], 'FontSize', fs);

subplot(2,4,5)
imshow(I_GL1);
title(['GL1 (psnr = ', num2str(psnr_GL1) ')'], 'FontSize', fs);
xlabel(['Time = ', num2str(time_GL1)], 'FontSize', fs);

subplot(2,4,6)
imshow(I_FGL1);
title(['FGL1 (psnr = ', num2str(psnr_FGL1) ')'], 'FontSize', fs);
xlabel(['Time = ', num2str(time_FGL1)], 'FontSize', fs);

subplot(2,4,7)
imshow(I_GLQ);
title(['GLQ (psnr = ', num2str(psnr_GLQ) ')'], 'FontSize', fs);
xlabel(['Time = ', num2str(time_GLQ)], 'FontSize', fs);

subplot(2,4,8)
imshow(I_FGLQ);
title(['FGLQ (psnr = ', num2str(psnr_FGLQ) ')'], 'FontSize', fs);
xlabel(['Time = ', num2str(time_FGLQ)], 'FontSize', fs);

%%
figure
imshow(img0); 
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_clean.pdf';

figure
imshow(I_noisy);
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_noisy.pdf';

figure
imshow(I_OMP);
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_OMP.pdf';

figure
imshow(I_IRLS);
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_IRLS.pdf';

figure
imshow(I_GL1);
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_GL1.pdf';

figure
imshow(I_FGL1);
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_FGL1.pdf';

figure
imshow(I_GLQ);
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_GLQ.pdf';

figure
imshow(I_FGLQ);
set(gcf, 'Color', 'w');
% export_fig 'figs_img_new/MRI11_k12_FGLQ.pdf';






%% Show Results

% figure(1)
% imshow(img0);
% export_fig 'figs_img/flower_GT.pdf'; % GT: Ground Truth
% 
% figure(2)
% imshow(X_OMP);
% export_fig 'figs_img/flower_OMP.pdf';
% 
% figure(3)
% imshow(X_IRLS);
% export_fig 'figs_img/flower_IRLS.pdf';
% 
% figure(4)
% imshow(X_GL2);
% export_fig 'figs_img/flower_GL2.pdf';
% 
% figure(5)
% imshow(X_GLQ);
% export_fig 'figs_img/flower_GLQ.pdf';
% 
% figure(6)
% imshow(X_FGLQ);
% export_fig 'figs_img/flower_FGLQ.pdf';


