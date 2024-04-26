function excution_blur()
image = imread('vandy .png');
result = blur(image,3);
imshow(result)
end