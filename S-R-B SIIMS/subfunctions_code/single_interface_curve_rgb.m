function I256 = single_interface_curve_rgb(mask_l1, mask_l2, T_or_empty, out_path, lineWidth)


    if nargin < 3 || isempty(T_or_empty), T_or_empty = []; end
    if nargin < 4, out_path = ''; end
    if nargin < 5, lineWidth = 1; end
    mask_l1 = logical(mask_l1);
    mask_l2 = logical(mask_l2);
    se   = strel('square',3);
    touch = imdilate(mask_l1,se) & imdilate(mask_l2,se);
    per1 = bwperim(mask_l1);
    per2 = bwperim(mask_l2);
    band = (per1 | per2) & touch;
    curve = bwmorph(per2,'thin',Inf);
    curve = imresize(curve,[256 256],'nearest');
    if lineWidth > 1
        se2 = strel('disk', floor(lineWidth/2), 0);
        curve = imdilate(curve, se2);
    end


    if isempty(T_or_empty)
        I256 = zeros(256,256,3);
    else
        if ndims(T_or_empty)==2
            bg = mat2gray(imresize(T_or_empty,[256 256]));
            I256 = repmat(bg,[1 1 3]);  
        else
            I256 = im2double(imresize(T_or_empty,[256 256]));
        end
    end


    R = I256(:,:,1); G = I256(:,:,2); B = I256(:,:,3);
    R(curve) = 1;
    G(curve) = 0;
    B(curve) = 0;
    I256 = cat(3,R,G,B);

end
