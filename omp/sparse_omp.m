% *******************************************************************************************
% Compressed Sensing Optimization Methods
% Alankar Kotwal <alankar.kotwal@iitb.ac.in>
% Sparse measurement matrix driver file
% Run the code after CDing to this file's directory.
% *******************************************************************************************

% Load the image
im = imread('../data/barbara256.png');
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
m = ceil(64*0.05);

% Get the DCT matrix
D = kron(dctmtx(8)', dctmtx(8)');

% Measurement matrix
randNums = randperm(64, m);
phim = zeros(m, 64);
for i=1:m
    phim(i, :) = D(randNums(i), :);
end
A = phim * D;

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

% OMP machao
op = zeros(imSize(1), imSize(2));
opMask = zeros(imSize(1), imSize(2));
sqPatchErr = 0;
for i=1:measSize(1)
    for j=1:measSize(2)
        
        disp([i j]);
        
        res = meas(:, i, j);
        suppSet = zeros(1, 64);
        suppSetIdx = 1;
        while(sum(res .* res) > 10000/m^2)
            
            coeffs = abs(sum(A .* repmat(res, 1, 64))./sum(A .* A));
            [~, idx] = max(coeffs);
            
            if(isempty(find(suppSet == idx, 1)))
                suppSet(suppSetIdx) = idx;
                suppSetIdx = suppSetIdx + 1;
            end
  
            Ai = zeros(m, suppSetIdx-1);
            for k=1:suppSetIdx-1
                Ai(:, k) = A(:, suppSet(k));
            end
            
            si = Ai\meas(:, i, j);
            
            res = meas(:, i, j);
            for k=1:suppSetIdx-1
                res = res - si(k)*Ai(:, k);
            end
            
        end

        patch = zeros(64, 1);
        for k=1:suppSetIdx-1
            patch = patch + si(k)*D(:, suppSet(k));
        end
        op(i:i+7, j:j+7) = op(i:i+7, j:j+7) + reshape(patch, 8, 8);
        opMask(i:i+7, j:j+7) = opMask(i:i+7, j:j+7) + ones(8);
        
        sqPatchErr = sqPatchErr + sum(sum((double(im(i:i+7, j:j+7)) - reshape(patch, 8, 8)).^2));
        
    end
end

output = op./opMask;
imshow(output/max(max(output)));

MSPE = sqPatchErr/(measSize(1)*measSize(2));
MSIE = sum(sum((output-double(im)).^2))/(imSize(1)*imSize(2));
disp([MSPE MSIE]);