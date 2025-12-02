function u = im2uint8(img)
    if islogical(img)
        u = uint8(img) * 255;
    elseif isa(img, 'uint8')
        u = img;
    elseif isa(img, 'uint16') || isa(img, 'double') || isa(img, 'single')
        u = uint8(min(max(round(img * 255), 0), 255));
    else
        error('Unsupported image type');
    end
end
