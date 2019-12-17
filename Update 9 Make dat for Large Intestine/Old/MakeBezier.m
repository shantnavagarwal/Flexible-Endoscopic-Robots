image = imread('Large Intestine sketch.png');
image = imresize(image, [320, 240]);
imshow(image)
figure()
bw = im2bw(image);
bw = imcomplement(bw);
[rows, col] = find(bw);
scatter(flip(rows), flip(col))
Q = [rows, col];
