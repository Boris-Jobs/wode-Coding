clc
clear

addpath('../');

img = double(imread('a.jpg'))/255;
img2=double(imread('a_smd.png'))/255;
img3 = double(imread('a.jpg'))/255;
ws = 8; 
a1=floor(size(img,1)/ws);
b1=floor(size(img,2)/ws);
img = img(1:a1*ws, 1:b1*ws);

no_patches1 = size(img, 1) / ws;
no_patches2 = size(img, 2) / ws;
X = zeros(2,ws^2);
k = 1;
for i = (1:no_patches1)
    for j = (1:no_patches2)
        r1 = (i-1)*ws+1:i*ws;
        r2 = (j-1)*ws+1:j*ws;
        patch = img(r1, r2);
        X(k,:) = patch(:);
        k = k + 1;
    end
end

% apply Robust PCA
lambda = 0.02; % close to the default one, but works better
tic
[L, S] = RobustPCA(X, lambda, 1.0, 1e-6);
toc

% reconstruct the image from the overlapping patches in matrix L
img_low_rank = zeros(size(img));
img_sparse = zeros(size(img));
k = 1;
for i = (1:no_patches1)
    for j = (1:no_patches2)
        patch = reshape(L(k,:), ws, ws);
        r1 = (i-1)*ws+1:i*ws;
        r2 = (j-1)*ws+1:j*ws;
        img_low_rank(r1, r2) = img_low_rank(r1, r2) + patch;
        patch = reshape(S(k,:), ws, ws);
        img_sparse(r1, r2) = img_sparse(r1, r2) + patch;
        k = k + 1;
    end
end

norm1=zeros(size(S,1),1);
for i=1:size(S,1)
    for j=1:size(S,2)
        norm1(i,1)=norm1(i,1)+abs(S(i,j));
    end
end

S1=zeros(size(S));
for i=1:size(S,1)
    if norm1(i)>=quantile(norm1,0.8,1)
        for j=1:size(S,2)
            S1(i,j)=255;
        end
    elseif (norm1(i)>=quantile(norm1,0.7,1))&&(norm1(i)<quantile(norm1,0.8,1))
        for j=1:size(S,2)
            S1(i,j)=175;
        end
    elseif(norm(i)>=quantile(norm1,0.6,1))&&(norm1(i)<quantile(norm1,0.7,1))
        for j=1:size(S,2)
            S1(i,j)=85;
        end
    else
        for j=1:size(S,2)
            S1(i,j)=0;
        end
    end
end

saliencemap=zeros(size(img));
k=1;
for i = (1:no_patches1)
    for j = (1:no_patches2)
        r1 = (i-1)*ws+1:i*ws;
        r2 = (j-1)*ws+1:j*ws;
        patch = reshape(S1(k,:), ws, ws);
        saliencemap(r1, r2) = saliencemap(r1, r2) + patch;
        k = k + 1;
    end
end

% show the results
figure;
subplot(2,3,1), imshow(img3,[]), title('Original')
subplot(2,3,2), imshow(img_low_rank,[]), title('low rank image')
bw_s = imbinarize(img_sparse,0.0000000001);
bw_s=bw_s*255;
subplot(2,3,3), imshow(bw_s), title('sparse image')
subplot(2,3,4),imshow(saliencemap,[]),title('salience map')
subplot(2,3,5),imshow(img2,[]),title('SMD result')
fprintf(1, 'ws=%d\tlambda=%f\trank(L)=%d\tcard(S)=%d\terr=%f\n', ...
    ws, lambda, rank(L), nnz(S), norm(img - img_low_rank, 'fro'));
