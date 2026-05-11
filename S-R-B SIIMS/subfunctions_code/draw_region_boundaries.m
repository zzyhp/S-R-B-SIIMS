function I256 = draw_region_boundaries(mask_l1, mask_l2, mask_l3, T_or_empty, out_path, lineWidth)

    if nargin < 4 || isempty(T_or_empty), T_or_empty = []; end
    if nargin < 5, out_path = ''; end
    if nargin < 6, lineWidth = 1; end
    lineWidth = max(1, min(10, round(lineWidth)));

    minArea = 5; 
    mask_l1 = bwareaopen(logical(mask_l1), minArea);
    mask_l2 = bwareaopen(logical(mask_l2), minArea);
    mask_l3 = bwareaopen(logical(mask_l3), minArea);

    se = strel('disk', 3);
    
    touch12 = imdilate(mask_l1, se) & imdilate(mask_l2, se);
    touch13 = imdilate(mask_l1, se) & imdilate(mask_l3, se);
    touch23 = imdilate(mask_l2, se) & imdilate(mask_l3, se);
    
    curve12 = bwperim(mask_l1, 8) & touch12; 
    curve13 = bwperim(mask_l1, 8) & touch13;
    curve23 = bwperim(mask_l2, 8) & touch23;

    minBoundaryArea = 3; 
    curve12 = bwareaopen(curve12, minBoundaryArea);
    curve13 = bwareaopen(curve13, minBoundaryArea);
    curve23 = bwareaopen(curve23, minBoundaryArea);

    if lineWidth > 1
        if lineWidth <= 5
            se_line = strel('disk', lineWidth-1, 0);
        else
            se_line = strel('rectangle', [lineWidth lineWidth]);
        end
        curve12 = imdilate(curve12, se_line);
        curve13 = imdilate(curve13, se_line);
        curve23 = imdilate(curve23, se_line);
    end

    targetSize = [256 256];
    if ~isequal(size(curve12), targetSize)
        curve12 = imresize(curve12, targetSize, 'nearest');
        curve13 = imresize(curve13, targetSize, 'nearest');
        curve23 = imresize(curve23, targetSize, 'nearest');
    end

    if isempty(T_or_empty)
        I256 = zeros(256, 256, 3, 'double');
    else
        if ndims(T_or_empty) == 2
            bg = mat2gray(imresize(T_or_empty, targetSize));
            I256 = repmat(bg, [1 1 3]);
        else
            I256 = im2double(imresize(T_or_empty, targetSize));
        end
    end

    R = I256(:,:,1); G = I256(:,:,2); B = I256(:,:,3);
    
    R(curve12) = 1; G(curve12) = 0; B(curve12) = 0; 
    R(curve13) = 0; G(curve13) = 1; B(curve13) = 0;
    R(curve23) = 0; G(curve23) = 0; B(curve23) = 1; 
    
    I256 = cat(3, R, G, B);
