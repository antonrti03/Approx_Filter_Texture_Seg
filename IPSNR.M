function [mse, rmse, psnr] = ipsnr(grayImage, noisyImage, BD)
%[mse, rmse, psnr] = ipsnr(grayImage, noisyImage)
% MSE: = rmse^2: Loi binh phuong trung binh
% f1: anh sau giai nen
% f2: anh truoc khi nen
%Tinh toan va hien thi loi giua anh truoc va sau khi giai nen
%Neu bo qua scale thi df=1.
%Compute the root-mean-square error.
%--------------------------------------------------------------------------
% e = double(f1) - double(f2);
% fd = double(f2);
% PSNR = 10*log10(sum(fd(:).^2)/sum(e(:).^2));
%--------------------------------------------------------------------------
% Now we have our two images and we can calculate the PSNR.
% First, calculate the "square error" image.
% Make sure they're cast to floating point so that we can get negative differences.
% Otherwise two uint8's that should subtract to give a negative number
% would get clipped to zero and not be negative.

squaredErrorImage = (double(grayImage) - double(noisyImage)).^ 2;

% Display the squared error image.
%subplot(2, 2, 3);
%imshow(squaredErrorImage, []);
%title('Squared Error Image', 'FontSize', fontSize);
% Sum the Squared Image and divide by the number of elements
% to get the Mean Squared Error.  It will be a scalar (a single number).

[rows, columns] = size(grayImage);

mse = sum(squaredErrorImage(:)) / (rows * columns);

rmse = sqrt(sum(squaredErrorImage(:)) / (rows * columns));

% Calculate PSNR (Peak Signal to Noise Ratio) from the MSE according to the formula.

psnr = 10 * log10( (2^BD-1)^2 / mse);

% Alert user of the answer.
%message = sprintf('The mean square error is %.2f.\nThe PSNR = %.2f', mse, PSNR);
%msgbox(message);
end
%%
%  function PSNR = ipsnr(distImg, origImg)
% 
%             origImg = double(origImg);
%             distImg = double(distImg);
% 
%             [M N] = size(origImg);
%             error = origImg - distImg;
%             MSE = sum(sum(error .* error)) / (M * N);
% 
%             if(MSE > 0)
%             PSNR = 10*log(255*255/MSE) / log(10);
%             else
%             PSNR = 99;
%             end 
% 
%  end
%%
%There is some mistakes with above code where it should consider log10 in Matlab. 
%So the correct code (even for color image) is as below:
% function PSNR = MyPsnr(distImg, origImg)
%         origImg = double(origImg);
%         distImg = double(distImg);
% 
%         [M N P] = size(origImg);
%         error = origImg - distImg;
%         MSE = sum(sum(sum(error.^2))) / (M * N * P);
% 
%         if(MSE > 0)
%             PSNR = 20*log10(max(max(max(origImg))))-10*log10(MSE);
%         else
%             PSNR = 99;
%         end 
% end
%%