% *******************************************************************************************
% Compressed Sensing Optimization Methods
% Alankar Kotwal <alankar.kotwal@iitb.ac.in>
% Random measurement matrix driver file
% Run the code after CDing to this file's directory.
% *******************************************************************************************

K = 10;
t = 1;
eps = 1000;
maxIter = 1000;

% Load the image
im = imread('../data/barbara256.png');
im = imresize(im, 0.25);
imSize = size(im);
measSize = imSize - 7;

% Extract patches from the Barbara image
patches = zeros(64, measSize(1), measSize(2));
for i=1:measSize(1)
    for j=1:measSize(2)
        
        patches(:, i, j) = reshape(im(i:i+7, j:j+7), 64, 1);
        
    end
end

% Parameters
m = ceil(64*0.4);

% Measurement matrix
phi = randn(64);
phim = phi(1:m, :);

% Extract measurements using this matrix
meas = zeros(m, measSize(1), measSize(2));
for i=1:measSize(1)
    for j=1:measSize(2)
        
        meas(:, i, j) = phim*patches(:, i, j);
        
    end
end

% Add some noise
sigma = 0.05*mean(mean(mean(abs(meas))));
meas = meas + sigma*randn(m, measSize(1), measSize(2));

% Get the DCT matrix
D = kron(dctmtx(8)', dctmtx(8)');
A = phim * D;

% SCA machao
op = zeros(imSize(1), imSize(2));
opMask = zeros(imSize(1), imSize(2));
sqPatchErr = 0;
for i=1:measSize(1)
    for j=1:measSize(2)
        
        disp([i j]);
        
        % Solve
        % A*z = meas(:, i, j)
        % such that you minimize l0 norm of z.
        z = solveSCA(A, meas(:, i, j), K, t, eps, maxIter);
        
        patch = D*z;
        
        op(i:i+7, j:j+7) = op(i:i+7, j:j+7) + reshape(patch, 8, 8);
        opMask(i:i+7, j:j+7) = opMask(i:i+7, j:j+7) + ones(8);
        
%         opt = op./opMask;
%         opt(isnan(opt)) = 0;
%         imshow(opt/max(max(opt)));
%         drawnow;
        
        sqPatchErr = sqPatchErr + sum(sum((double(im(i:i+7, j:j+7)) - reshape(patch, 8, 8)).^2));
        
    end
end

output = op./opMask;
%imshow(output/max(max(output)));
imwrite('output-'+int2str(now));

MSPE = sqPatchErr/(measSize(1)*measSize(2));
MSIE = sum(sum((output-double(im)).^2))/(imSize(1)*imSize(2));
disp([MSPE MSIE]);
