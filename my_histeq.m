function [output_img] = my_histeq(img)

    [rows, cols] = size(img);
    num_pixels = rows * cols;
    hist_counts = my_imhist(img);
    
    pdf = hist_counts / num_pixels;
    
    cdf = cumsum(pdf);
    
    equalized_values = round(cdf * 255);
    
    output_img = uint16process(zeros(size(img)));
    for i = 1:rows
        for j = 1:cols
            output_img(i,j) = equalized_values(img(i,j)+1);
        end
    end
end