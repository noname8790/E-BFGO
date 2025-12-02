function gray_img = custom_mat2gray(img)
    min_val = min(img(:));
    max_val = max(img(:));
    gray_img = (double(img) - min_val) / (max_val - min_val);
end