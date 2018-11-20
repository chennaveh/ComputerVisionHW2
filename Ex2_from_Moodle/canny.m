function [output] = canny(file_name, sigma ,L_th, H_th)
% Canny Edge Detector

% calculate derivatives
mask_size = round(2 * sigma) + 1;
G_dx = Deriv_Gauss_x(sigma, mask_size);
G_dy = Deriv_Gauss_y(sigma, mask_size);

% generate I_x, I_y
I = imread(file_name);
I_x = conv2(I, G_dx,'same');
I_y = conv2(I, G_dy,'same');

% generate G_orientation and G_magnitute
G_magnitute = (I_x.^2 + I_y.^2).^0.5;
G_orientation = atan(I_y ./ I_x);

% generate Et for non-maximum suppression
Et = thinning(G_magnitute, G_orientation);

% hysteresis thresholding - 8 neighbors
im_uint8 = uint8(Et / max(Et(:)) * 255);
bin_im_lo_th = (im_uint8 >= L_th) & (im_uint8 < H_th); % mid values 
bin_im_hi_th = im_uint8 >= H_th;                       % high values 
mask = true(3, 3);
bin_hi_neighbors = logical(conv2(bin_im_hi_th, mask)); % find neighbores of high values
bin_result = (bin_hi_neighbors(2:end-1, 2:end-1) & bin_im_lo_th) | bin_im_hi_th; % gather requested pixels

output = double(bin_result) .* Et;

end

function [output] = Deriv_Gauss_x(sigma, mask_size)
% generate one directional gaussian derivative
x = meshgrid(linspace(-((mask_size)/2 ), (mask_size)/2 , mask_size));
output = -x.*(exp(-(x.^2)/(2*sigma^2))).*(sigma^4);
% output = output./sum(output(:));
end

function [output] = Deriv_Gauss_y(sigma, mask_size)
% generate one derictional gaussian window using transpose
output = Deriv_Gauss_x(sigma, mask_size)';
end

function [output] = thinning(G_strength, G_orientation)
% round orientation to [0, pi/4, pi/2, 3*pi/4]
oriented = mod(G_orientation + pi, pi);
oriented = mod(round(oriented / pi * 4), 4);
oriented = oriented * pi / 4;

% generate new image with enhanced values
size_new_im = size(G_strength) + 2;
padded = zeros(size_new_im); % pad zero frame
padded(2:end-1, 2:end-1) = G_strength;
output = zeros(size(G_strength));
for i = 2 : size_new_im(1)-1
    for j = 2 : size_new_im(2)-1
        if     (oriented(i-1, j-1) == 0      && padded(i,j) >= padded(i,j+1)   && padded(i,j) >= padded(i,j-1))
            output(i-1,j-1) = G_strength(i-1,j-1);
        elseif (oriented(i-1, j-1) == pi/4   && padded(i,j) >= padded(i-1,j-1) && padded(i,j) >= padded(i+1,j+1))
            output(i-1,j-1) = G_strength(i-1,j-1);
        elseif (oriented(i-1, j-1) == pi/2   && padded(i,j) >= padded(i+1,j)   && padded(i,j) >= padded(i-1,j))
            output(i-1,j-1) = G_strength(i-1,j-1);
        elseif (oriented(i-1, j-1) == 3*pi/4 && padded(i,j) >= padded(i-1,j+1) && padded(i,j) >= padded(i+1,j-1))
            output(i-1,j-1) = G_strength(i-1,j-1);
        end
    end
end

end







