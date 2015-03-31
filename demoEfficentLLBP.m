inImg=imread('peppers.png');
%inImg = rgb2gray(inImg);
filtUpLLBP = generateFilterLLBP(5);
nFiltDims = [7, 7];
hLBP = @efficientLBP;       % @pixelwiseLBP, @efficientLBP
isScale = true;

tic; 
filtLLBP = lineFiltLBP(inImg, nFiltDims, hLBP, isScale);
timeFiltLBP = toc;

tic;
shiftLLBP = lineShiftLBP(inImg, nFiltDims, isScale);
timeShiftLBP = toc;

figure;
subplot(1, 3, 1);
imshow(inImg);  
title('Input image'); 

subplot(1, 3, 2); 
imshow(filtLLBP);  
title( sprintf('Filter based LLBP, runtime %.1f [sec]', timeFiltLBP) );

subplot(1, 3, 3); 
imshow(shiftLLBP);  
title( sprintf('Filter based LLBP, runtime %.1f [sec]', timeShiftLBP) );
