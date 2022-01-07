%% Please run the code section by section
%% 1. Read a jpg file and convert it to a grayscale image to save more storage space
clear all, close all, clc
A = imread('Jin_Hoontae.jpg');
B = rgb2gray(A);

figure();
subplot(1,2,1)
imshow(A) %Original Image
title('Original image')

subplot(1,2,2)
imshow(B) %Grayscale Image
title('Gray image')
%% 2. Perform FFT without inbuilt functions
% 2.1. Discrete Fourier Transform Construction for 2D Matrix
[Row,Col] = size(B);
Row2 = zeros(Row, Row);
Col2 = zeros(Col, Col);

for n = 0 : (Row - 1)
    for k = 0 : (Row - 1)
        Row2(n+1, k+1) = exp((2 * pi * 1i) * (k * n / Row));
    end    
end

for n2 = 0 : (Col - 1)
    for k2 = 0 : (Col - 1)
        Col2(k2+1, n2+1) = exp((2 * pi * 1i) * (k2 * n2 / Col));
    end    
end

Bt = Row2*double(B)*Col2; %Discrete Fourier Transform Equation for 2D matrix

% 2.2. FFT and FFT Shift Construction to obtain the spectrum of centered frequency coefficients
sz = ceil(size(Bt)/2); 
%Line 38-39 shifts the low frequency coefficients located at each edge to the center for the better visualization
Bt2 = Bt([sz(1)+1:end, 1:sz(1)], [sz(2)+1:end, 1:sz(2)]);
Blog = log(abs(Bt2)+1); %Computation on Log-scale for better analysis - FFT

figure();
subplot(1,2,1)
imshow(mat2gray(Blog),[]) %Shows the computated frequency domain by FFT
title('without built-in functions')
subplot(1,2,2)
graylog = log(abs(fftshift(fft2(B)))+1);
imshow(mat2gray(graylog),[]); % frequency spectrum with built-in functions
title('with built-in functions')
%% 3.Image Compression (Inverse DFT and Eliminiation of small coefficients)
% 3.1. Inverse Fourier Transform Construction for 2D Matrix
[Row, Col] = size(Bt);
Row2        = zeros(Row, Row);
Col2        = zeros(Col, Col);
 
for n = 0 : (Row - 1)
    for k = 0 : (Row - 1)
        Row2(n+1,k+1) = exp(-2 * pi * 1i / Row * n * k);
    end    
end
 
for n2 = 0 : (Col - 1)
    for k2 = 0 : (Col - 1)
        Col2( k2+1,n2+1) = exp(-2 * pi * 1i / Col * n2 * k2);
    end    
end
% Bi= (Row2*double(Bt)*Col2)/(Row*Col) : Inverse Discrete Fourier Transform Equation.
% This will be used after zeroing out small coefficients in the frequency domain 
% to produce compressed images.

% 3.2. Image compression by 30%, 5%, 1% and 0.2% for comparison
Bsort = sort(abs(Bt(:))); % Sort by magnitude
count = 1;
figure();
for keep=[.3 .05 .01 .002]
    subplot(2,2,count)
    thresh = Bsort(floor((1-keep)*length(Bsort)));
    ind = abs(Bt) > thresh; %Zero out small coefficients
    NumofCoef = Row*Col - sum(sum(ind)); %Count how many coefficients are removed
    percent = 100 - NumofCoef/(Row*Col)*100;
    Bt_new = Bt.*ind;
    Bi = (Row2*double(Bt_new)*Col2)/(Row*Col); %Inverse Discrete Fourier Transform with small coefficients emliminated
    Comp_Bi = uint8(round(abs(Bi))); %FFT-Compressed Images
    imshow(Comp_Bi)
    title(['',num2str(percent),'%', ' Compression Ratio'],'Fontsize',10)
    count = count + 1;
end
%% 4. FFT Frequency Coefficients comparison
Bsort = sort(abs(Bt(:)));
counter = 1;
figure();
for keep=[.3 .05 .01 .002] 
    subplot(2,2,counter)
    thresh = Bsort(floor((1-keep)*length(Bsort)));
    ind = abs(Bt) > thresh; 
    NumofCoef = Row*Col - sum(sum(ind));
    percent = 100 - NumofCoef/(Row*Col)*100;
    Bt_new = Bt.*ind;
    sz = ceil(size(Bt_new)/2); 
    Bt_Centered = Bt_new([sz(1)+1:end, 1:sz(1)], [sz(2)+1:end, 1:sz(2)]); %FFT Shift
    Blog = log(abs(Bt_Centered)+1); %FFT 
    imshow(mat2gray(Blog))
    title(['',num2str(percent),'%'],'Fontsize',10)
    counter = counter + 1;
end
%% 5. Pad the original image before applying Gaussian Low-Pass filter
[m n] = size(B);
P = 2*m;
Q = 2*n;
Padded_img = zeros(P,Q);
Padded_img(1:m,1:n) = B;
for i = 1:P
    for j = 1:Q
        Padded_img(i,j) = Padded_img(i,j)*(-1)^(i+j);
    end
end
%% 6. Perform FFT without inbuilt functions (Same as the step 2)
% 6.1. Discrete Fourier Transform Construction for 2D Matrix
[Row,Col] = size(Padded_img);
Row2 = zeros(Row, Row);
Col2 = zeros(Col, Col);

for n = 0 : (Row - 1)
    for k = 0 : (Row - 1)
        Row2(n+1, k+1) = exp((2 * pi * 1i) * (k * n / Row));
    end    
end

for n2 = 0 : (Col - 1)
    for k2 = 0 : (Col - 1)
        Col2(k2+1, n2+1) = exp((2 * pi * 1i) * (k2 * n2 / Col));
    end    
end

Bt_Gaus = Row2*double(Padded_img)*Col2; %Discrete Fourier Transform Equation for 2D matrix
Blog_Gaus = log(abs(Bt_Gaus)+1); %Computation on Log-scale for better analysis - FFT
%% 7. Construct Guassian Low-Pass filter
[m n] = size(Blog_Gaus);
img_frame = zeros(m,n);
D = zeros(m,n);

for i = 1:m
    for j = 1:n
        D(i,j) = [(i-(P/2))^2 + (j-(Q/2))^2]^(1/2);
        img_frame(i,j) = exp(-D(i,j)^2/(2*180^2)); %Where 180 is the cutoff frequency value
    end
end
figure();
imshow(mat2gray(img_frame))
title('Gaussian Low-Pass Filter')
%% 8. Convolute the frequency spectrum with the Gaussian filter (to smooth the image)
G = img_frame.*Bt_Gaus;
[Row, Col] = size(G);
Row2        = zeros(Row, Row);
Col2        = zeros(Col, Col);
 
for n = 0 : (Row - 1)
    for k = 0 : (Row - 1)
        Row2(n+1,k+1) = exp(-2 * pi * 1i / Row * n * k);
    end    
end
 
for n2 = 0 : (Col - 1)
    for k2 = 0 : (Col - 1)
        Col2(k2+1,n2+1) = exp(-2 * pi * 1i / Col * n2 * k2);
    end    
end
Bi_Gaus = (Row2*G*Col2)/(Row*Col);

for i = 1:P
    for j = 1:Q
        Bi_Gaus(i,j) = Bi_Gaus(i,j)*(-1)^(i+j); % Shift the centered low frequencies back to the original locations
    end
end
Comp_Bi_Gaus = (real((Bi_Gaus)));
Comp_Bi_Gaus = Comp_Bi_Gaus(1:size(B,1),1:size(B,2)); % Remove the padded area

figure();
subplot(1,2,1)
imshow(mat2gray(Comp_Bi_Gaus))
title('Gaussian-filtered image')

subplot(1,2,2)
imshow(B)
title('Original Image')
%% 9. Re-perform FFT on the Gaussian-compressed image (same as the step 2)
[Row,Col] = size(Comp_Bi_Gaus);
Row2 = zeros(Row, Row);
Col2 = zeros(Col, Col);

for n = 0 : (Row - 1)
    for k = 0 : (Row - 1)
        Row2(n+1, k+1) = exp((2 * pi * 1i) * (k * n / Row));
    end    
end

for n2 = 0 : (Col - 1)
    for k2 = 0 : (Col - 1)
        Col2(k2+1, n2+1) = exp((2 * pi * 1i) * (k2 * n2 / Col));
    end    
end

Bt_Gaus = Row2*double(Comp_Bi_Gaus)*Col2;

sz_Gaus = ceil(size(Bt_Gaus)/2); 
Bt2_Gaus = Bt_Gaus([sz_Gaus(1)+1:end, 1:sz_Gaus(1)], [sz_Gaus(2)+1:end, 1:sz_Gaus(2)]);
Blog_Gaus = log(abs(Bt2_Gaus)+1);

figure();
subplot(1,2,1)
imshow(mat2gray(Blog_Gaus),[]) 
title('without built-in functions')
subplot(1,2,2)
graylog = log(abs(fftshift(fft2(Comp_Bi_Gaus)))+1);
imshow(mat2gray(graylog),[]);
title('with built-in functions')
%% 10.Image Compression (Same as the step 3)
[Row, Col] = size(Bt_Gaus);
Row2        = zeros(Row, Row);
Col2        = zeros(Col, Col);
 
for n = 0 : (Row - 1)
    for k = 0 : (Row - 1)
        Row2(n+1,k+1) = exp(-2 * pi * 1i / Row * n * k);
    end    
end
 
for n2 = 0 : (Col - 1)
    for k2 = 0 : (Col - 1)
        Col2( k2+1,n2+1) = exp(-2 * pi * 1i / Col * n2 * k2);
    end    
end

Bsort_Gaus = sort(abs(Bt_Gaus(:)));
count = 1;
figure();
for keep=[.3 .05 .01 .002]
    subplot(2,2,count)
    thresh_Gaus = Bsort_Gaus(floor((1-keep)*length(Bsort_Gaus)));
    ind_Gaus = abs(Bt_Gaus) > thresh_Gaus; 
    NumofCoef_Gaus = Row*Col - sum(sum(ind_Gaus)); 
    percent_Gaus = 100 - NumofCoef_Gaus/(Row*Col)*100;
    Bt_new_Gaus = Bt_Gaus.*ind_Gaus;
    Bi_Gaus = (Row2*double(Bt_new_Gaus)*Col2)/(Row*Col); 
    Comp_Bi_Gaus = uint8(round(abs(Bi_Gaus)));
    imshow(Comp_Bi_Gaus)
    title(['',num2str(percent_Gaus),'%', ' Compression Ratio'],'Fontsize',10)
    count = count + 1;
end
%% 11. Frequency Coefficients comparison (Same step as the step 4)
Bsort_Gaus = sort(abs(Bt_Gaus(:)));
counter = 1;
figure();
for keep=[.3 .05 .01 .002] 
    subplot(2,2,counter)
    thresh_Gaus = Bsort_Gaus(floor((1-keep)*length(Bsort_Gaus)));
    ind_Gaus = abs(Bt_Gaus) > thresh_Gaus; 
    NumofCoef_Gaus = Row*Col - sum(sum(ind_Gaus));
    percent_Gaus = 100 - NumofCoef_Gaus/(Row*Col)*100;
    Bt_new_Gaus = Bt_Gaus.*ind_Gaus;
    sz = ceil(size(Bt_new_Gaus)/2); 
    Bt_Centered_Gaus = Bt_new_Gaus([sz(1)+1:end, 1:sz(1)], [sz(2)+1:end, 1:sz(2)]); %FFT Shift
    Blog_Gaus = log(abs(Bt_Centered_Gaus)+1); %FFT 
    imshow(mat2gray(Blog_Gaus))
    title(['',num2str(percent_Gaus),'%'],'Fontsize',10)
    counter = counter + 1;
end
%% 12. Perform SVD without inbuilt functions
B = double(rgb2gray(A)); %Input Image
[Row,Col] = size(B);

% SVD Matrix Construction (B = U*S*V')
if  Row > Col
  
    [V_Evec,V_Eval] = eig(B'*B);
   
    for i = 1:Col
        V_norm(i) = norm(V_Evec(:,i));
        V(:,i) = (V_Evec(:,i))/V_norm(i); % Normalized Eigenvectors of V
    end
    
    S = zeros(size(B));
    s = fliplr(sqrt(V_Eval));
    s = flip(s,1); %Set the values in descending order.
    S(1:Col,:) = S(1:Col,:) + s(1:Col,:); %Sigma Matrix Construction
    
    [~,I] = sort(diag(V_Eval),'descend');
    V = V(:,I); %V Construction
    U = rdivide(B*V,diag(s)'); %U Construction
    
else %If Row <= Col
    
    [U_Evec,U_Eval] = eig(B*B');
    
    for i = 1:Row
        U_norm(i) = norm(U_Evec(:,i));
        U(:,i) = (U_Evec(:,i))/U_norm(i);
    end
    
    S = fliplr(sqrt(U_Eval));
    S = flip(S,1); %Use fliplr and flip functions to set the matrix in descending order
    
    [~,I] = sort(abs(diag(U_Eval)),'descend');
    U = U(:,I);    
    V = rdivide(B'*U,diag(S)');
end

count = 1;
figure();
for r = [10 15 20 25 30 35 40 45 50]; % r: Low rank value 
    Br = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    subplot(3,4,count), count = count + 1;
    imshow(mat2gray(Br)), axis off;
    title(['r=',num2str(r,'%d'), ', ',num2str(100*r*(Row+Col)/(Row*Col),'%2.2f'),'% storage'])
    
end
%% 13. Distortion Computation for FFT-compressed and SVD-compressed images
figure();
% The formula: D = norm(B-O)^2/norm(B)^2 where B is the original gray image
% and O is the compressed image.

% 13.1. Distortion Error of FFT-compressed image
[Row Col] = size(B);
for keep=0.002:0.02:0.3
    thresh = Bsort(floor((1-keep)*length(Bsort)));
    ind = abs(Bt) > thresh; 
    NumofCoef = Row*Col - sum(sum(ind));
    Bt_new = Bt.*ind;
    percent = 100 - NumofCoef/(Row*Col)*100;
    Bi = (Row2*double(Bt_new)*Col2)/(Row*Col); 
    Comp_Bi = (round(abs(Bi)));
    C = double(rgb2gray(A));
    FFT_Distortion = norm(C-Comp_Bi,'fro').^2/norm(C,'fro').^2;
    plot_FFT=plot(percent/100, 100*FFT_Distortion, 'k.'); hold on;
end

% 13.2. Distortion Error of FFT-compressed+Gaussian filtered image
for keep=0.002:0.02:0.3
    thresh_Gaus = Bsort_Gaus(floor((1-keep)*length(Bsort_Gaus)));
    ind_Gaus = abs(Bt_Gaus) > thresh_Gaus; 
    NumofCoef_Gaus = Row*Col - sum(sum(ind_Gaus));
    Bt_new_Gaus = Bt_Gaus.*ind_Gaus;
    percent_Gaus = 100 - NumofCoef_Gaus/(Row*Col)*100;
    Bi_Gaus = (Row2*double(Bt_new_Gaus)*Col2)/(Row*Col); 
    Comp_Bi_Gaus = round(abs(Bi_Gaus));
    C = double(rgb2gray(A));
    FFT_Distortion_Gaus = norm(C-Comp_Bi_Gaus,'fro').^2/norm(C,'fro').^2;
    plot_FFT_G=plot(percent_Gaus/100, 100*FFT_Distortion_Gaus, 'r.'); hold on;
end

% 13.3. Distortion Error of SVD-compressed images
for r = 5:5:100;
    Br = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    SVD_Comp_ratio = r*(Row+Col)/Row/Col;
    SVD_Distortion = norm(double(B)-Br,'fro').^2/norm(double(Br),'fro').^2;
    SVD_Distortion = 100*SVD_Distortion;
    plot_SVD=plot(SVD_Comp_ratio,SVD_Distortion,'b.'); hold on
end
legend([plot_FFT plot_FFT_G plot_SVD],{'FFT Compression','FFT Compression with Gaussian filter at 180 cutoff frequency','SVD Compression'})
title('Comparison between FFT and SVD for image compression')
xlabel('Compression ratio')
ylabel('Distortion percent')

%% 14. PSNR (Peak signal-to-noise ratio)
figure();
% 14.1. PSNR values for FFT-Compressed images
for keep=0.002:0.02:0.3
    thresh = Bsort(floor((1-keep)*length(Bsort)));
    ind = abs(Bt) > thresh; 
    NumofCoef = Row*Col - sum(sum(ind));
    Bt_new = Bt.*ind;
    percent = 100 - NumofCoef/(Row*Col)*100;
    Bi = (Row2*double(Bt_new)*Col2)/(Row*Col); 
    Comp_Bi = round(abs(Bi));
    
    % Mean Squared Error
    sz = size(B(:));
    I1 = B(:);
    I2 = Comp_Bi(:);
    FFT_MSE = sum((I1-I2).^2)/sz(1);
    % PSNR
    FFT_PSNR = 10*log10(255^2/FFT_MSE);
    plot_FFT=plot(percent/100,FFT_PSNR, 'k.'); hold on;
end

% 14.2. PSNR values for FFT-Compressed+Gaussian filtered images
for keep=0.002:0.02:0.3
    thresh_Gaus = Bsort_Gaus(floor((1-keep)*length(Bsort_Gaus)));
    ind_Gaus = abs(Bt_Gaus) > thresh_Gaus; 
    NumofCoef_Gaus = Row*Col - sum(sum(ind_Gaus));
    Bt_new_Gaus = Bt_Gaus.*ind_Gaus;
    percent_Gaus = 100 - NumofCoef_Gaus/(Row*Col)*100;
    Bi_Gaus = (Row2*double(Bt_new_Gaus)*Col2)/(Row*Col); 
    Comp_Bi_Gaus = round(abs(Bi_Gaus));
    
    % Mean Squared Error
    sz = size(B(:));
    I1 = B(:);
    I2 = Comp_Bi_Gaus(:);
    FFT_MSE_Gaus = sum((I1-I2).^2)/sz(1);
    % PSNR
    FFT_PSNR_Gaus = 10*log10(255^2/FFT_MSE_Gaus);
    plot_FFT_G=plot(percent_Gaus/100,FFT_PSNR_Gaus, 'r.'); hold on;
end

% 14.3. PSNR Values for SVD-Compressed images
for r = 5:5:100;
    By = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    SVD_Comp_ratio = r*(Row+Col)/Row/Col;
    SVD_PSNR = psnr(uint8(B),uint8(By));
    
    % Mean Squared Error
    sz = size(B(:));
    I1 = B(:);
    I2 = By(:);
    SVD_MSE = sum((I1-I2).^2)/sz(1);
    % PSNR
    SVD_PSNR = 10*log10(255^2/SVD_MSE);
    plot_SVD=plot(SVD_Comp_ratio,SVD_PSNR,'b.'); hold on
end
legend([plot_FFT plot_FFT_G plot_SVD],{'FFT Compression','FFT Compression with Gaussian filter at 180 cutoff frequency','SVD Compression'},'Location','southeast')
title('PNSR Analysis between FFT and SVD compressed images')
xlabel('Compression ratio')
ylabel('PSNR (dB)')