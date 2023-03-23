% Image Reconstruction Based on Compressive Sensing using Optimized Sensing Matrix
% Linear Algebra and Optimization 

% Linear Algebra Concepts - 
            % Sparse Matrix,
            % Matrix operations between Sparse Matrix and Measurement Matrix.

 % Optimization Concepts - 
            % Simulated Annealing for reducing the coherence between Sparse and Measurement Matrices - Future Work 
            % Orthogonal Matching Pursuit Algorithm for solving sparse signal and reconstructing the image.

% Clearing the previous workspaces, work variables and figures.            
tic
close all                                                                                                                                                                                                 
clear all
clc                                                                                                                                               

%Read Input
Original_Image = im2double(imread('Cars.jpg'));

% Converting the Original RGB Image into GrayScale Image - Since, RGB Image is Three Dimensional and Grayscale is Two Dimensional
Gray_Image=im2gray(Original_Image)                                                                                                                                                                                                          

% Generating a DCT Basis for calculating the Sparse Image
DCT_Basis = dctmtx(size(Gray_Image,1));

% Generating the Sparse Image using DCT Basis
Sparse_Image = DCT_Basis*Gray_Image*DCT_Basis';

% Displaying the Images
subplot(2,2,1)
imshow(Original_Image)
title('Original Image')

subplot(2,2,2);
imshow(Gray_Image)
title('Gray Scale Image')

subplot(2,2,3);
imshow(Sparse_Image);
title('Sparse Representation')

%Generate the Sparse Vector from the Sparse Representation
Gray_Size = size(Gray_Image);

% Random Number 
%N = 6300.1
N=1000;
% Reshape Sparse matrix into Vector by stacking columns under one another
Temp_Values=1;  
Sparse_Vector = zeros(Gray_Size(1,1)*Gray_Size(1,2),1);
for j = 1:Gray_Size(1,2)
    for i = 1:Gray_Size(1,1)
        Sparse_Vector(Temp_Values,1)=Sparse_Image(i,j);
        Temp_Values = Temp_Values+1;
    end
end

% Now, we are creating a Random Gaussian Sampling Matrix (i.e. Random Measurement Matrix)
Regeneration_Samples = 5000;
Random_Gaussian_Sampling_Matrix_Phi = rand(Regeneration_Samples,10*N);
for x = 1:Regeneration_Samples
    for y =1:10*N
        if (Random_Gaussian_Sampling_Matrix_Phi(x,y)+ 0.001 >= 1)
            Random_Gaussian_Sampling_Matrix_Phi(x,y) = 1;
        else
            Random_Gaussian_Sampling_Matrix_Phi(x,y) = 0;
        end
    end
end

% Random Measurement Process
Random_Measurement_Process = Random_Gaussian_Sampling_Matrix_Phi * Sparse_Vector; 

% Coherence Theta between Measurement Matrix and Sparse Basis

Theta = Random_Gaussian_Sampling_Matrix_Phi * dctmtx(size(10*N,1));

% Utilizing the Orthogonal Matching Pursuit (OMP) algorithm for reconstructing the gray scale image.

Final_Solution = (OMP(Random_Measurement_Process,Theta,10*N)');

% Converting Image vector back to a Pixel Matrix
% Defining an Empty Pixel Matrix
Solution_Pixel_Matrix = zeros(100,100); 
k=1;
j=1;
for i = 1 : 10000
    if(k > 100)
        k=1;
        j=j+1;
    end
    Solution_Pixel_Matrix(k,j) = Final_Solution(i);
    k=k+1;
end

% Utilizing the Invert 2D discrete Cosine Transform for generating image from the Pixel Matrix
Solution_Pixel_Matrix = idct2(Solution_Pixel_Matrix); 
subplot(2,2,4);
imshow(Solution_Pixel_Matrix)
title('Reconstructed Image')

% Calculating the Reconstruction Error by comparing the Original Gray Scale Image and Reconstructed Image
Reconstruction_Error = norm(Solution_Pixel_Matrix - Gray_Image) / norm(Gray_Image);

% Finding the Peak Signal-to-Noise Ratio for the reconstructed image

MSE=mean((Gray_Image(:)-Solution_Pixel_Matrix(:)).^2);
PSNR=10*log10(255^2/MSE);
disp(PSNR)
toc