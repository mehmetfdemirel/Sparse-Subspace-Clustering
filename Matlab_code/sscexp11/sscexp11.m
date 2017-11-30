%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      MAKE DATA                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 360 * 640;
Drows = dim;
N1_size = 41;

Y0mat = reshape(rgb2gray(im2double(imread('shift0.png'))), [dim,1]);
Y1mat = reshape(rgb2gray(im2double(imread('shift1.png'))), [dim,1]);
Y2mat = reshape(rgb2gray(im2double(imread('shift2.png'))), [dim,1]);
Y3mat = reshape(rgb2gray(im2double(imread('shift3.png'))), [dim,1]);
Y4mat = reshape(rgb2gray(im2double(imread('shift4.png'))), [dim,1]);
Y5mat = reshape(rgb2gray(im2double(imread('shift5.png'))), [dim,1]);
Y6mat = reshape(rgb2gray(im2double(imread('shift6.png'))), [dim,1]);
Y7mat = reshape(rgb2gray(im2double(imread('shift7.png'))), [dim,1]);
Y8mat = reshape(rgb2gray(im2double(imread('shift8.png'))), [dim,1]);
Y9mat = reshape(rgb2gray(im2double(imread('shift9.png'))), [dim,1]);
Y10mat = reshape(rgb2gray(im2double(imread('shift10.png'))), [dim,1]);
Y11mat = reshape(rgb2gray(im2double(imread('shift11.png'))), [dim,1]);
Y12mat = reshape(rgb2gray(im2double(imread('shift12.png'))), [dim,1]);
Y13mat = reshape(rgb2gray(im2double(imread('shift13.png'))), [dim,1]);
Y14mat = reshape(rgb2gray(im2double(imread('shift14.png'))), [dim,1]);
Y15mat = reshape(rgb2gray(im2double(imread('shift15.png'))), [dim,1]);
Y16mat = reshape(rgb2gray(im2double(imread('shift16.png'))), [dim,1]);
Y17mat = reshape(rgb2gray(im2double(imread('shift17.png'))), [dim,1]);
Y18mat = reshape(rgb2gray(im2double(imread('shift18.png'))), [dim,1]);
Y19mat = reshape(rgb2gray(im2double(imread('shift19.png'))), [dim,1]);
Y20mat = reshape(rgb2gray(im2double(imread('shift20.png'))), [dim,1]);
Y21mat = reshape(rgb2gray(im2double(imread('shift21.png'))), [dim,1]);
Y22mat = reshape(rgb2gray(im2double(imread('shift22.png'))), [dim,1]);
Y23mat = reshape(rgb2gray(im2double(imread('shift23.png'))), [dim,1]);
Y24mat = reshape(rgb2gray(im2double(imread('shift24.png'))), [dim,1]);
Y25mat = reshape(rgb2gray(im2double(imread('shift25.png'))), [dim,1]);
Y26mat = reshape(rgb2gray(im2double(imread('shift26.png'))), [dim,1]);
Y27mat = reshape(rgb2gray(im2double(imread('shift27.png'))), [dim,1]);
Y28mat = reshape(rgb2gray(im2double(imread('shift28.png'))), [dim,1]);
Y29mat = reshape(rgb2gray(im2double(imread('shift29.png'))), [dim,1]);
Y30mat = reshape(rgb2gray(im2double(imread('shift30.png'))), [dim,1]);
Y31mat = reshape(rgb2gray(im2double(imread('shift31.png'))), [dim,1]);
Y32mat = reshape(rgb2gray(im2double(imread('shift32.png'))), [dim,1]);
Y33mat = reshape(rgb2gray(im2double(imread('shift33.png'))), [dim,1]);
Y34mat = reshape(rgb2gray(im2double(imread('shift34.png'))), [dim,1]);
Y35mat = reshape(rgb2gray(im2double(imread('shift35.png'))), [dim,1]);
Y36mat = reshape(rgb2gray(im2double(imread('shift36.png'))), [dim,1]);
Y37mat = reshape(rgb2gray(im2double(imread('shift37.png'))), [dim,1]);
Y38mat = reshape(rgb2gray(im2double(imread('shift38.png'))), [dim,1]);
Y39mat = reshape(rgb2gray(im2double(imread('shift39.png'))), [dim,1]);
Y40mat = reshape(rgb2gray(im2double(imread('shift40.png'))), [dim,1]);


YY1mat = [Y0mat Y1mat Y2mat Y3mat Y4mat Y5mat Y6mat Y7mat Y8mat Y9mat Y10mat Y11mat Y12mat Y13mat Y14mat Y15mat Y16mat Y17mat Y18mat Y19mat Y20mat Y21mat Y22mat Y23mat Y24mat Y25mat Y26mat Y27mat Y28mat Y29mat Y30mat Y31mat Y32mat Y33mat Y34mat Y35mat Y36mat Y37mat Y38mat Y39mat Y40mat];

% Z1mat =  reshape([5,6,7,8,9;15,16,17,18,19;20,21,22,23,24],[15,1]);
% Z2mat =  reshape([6,7,8,9,15;16,17,18,19,20;21,22,23,24,5],[15,1]);
% Z3mat =  reshape([7,8,9,15,16;17,18,19,20,21;22,23,24,5,6],[15,1]);
% Z4mat =  reshape([8,9,15,16,17;18,19,20,21,22;23,24,5,6,7],[15,1]);
% Z5mat =  reshape([9,15,16,17,18;19,20,21,22,23;24,5,6,7,8],[15,1]);
% Z6mat =  reshape([15,16,17,18,19;20,21,22,23,24;5,6,7,8,9],[15,1]);
% Z7mat =  reshape([16,17,18,19,20;21,22,23,24,5;6,7,8,9,15],[15,1]);
% Z8mat =  reshape([17,18,19,20,21;22,23,24,5,6;7,8,9,15,16],[15,1]);
% Z9mat =  reshape([18,19,20,21,22;23,24,5,6,7;8,9,15,16,17],[15,1]);
% Z10mat = reshape([19,20,21,22,23;24,5,6,7,8;9,15,16,17,18],[15,1]);

Ymat = [YY1mat];



%%% Ymat = [Z1mat Z2mat Z3mat Z4mat Z5mat Z6mat Z7mat Z8mat Z9mat Z10mat];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    INITILIZATION                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = Drows;
N = N1_size;
maxIter = 10000;
C = zeros(N,N);
A = zeros(N,N);

E = zeros(D,N);

lowDelta = zeros(N,1);
delta = zeros(N,N);
rho = 1000;
lambdaZ = 0.001;
lambdaE = 0;

onevector = ones(N,1);
onematrix = onevector * transpose(onevector);
epsilon = 10 ^ (-4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    ITERATIVE ADMM                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=0:maxIter
    k
    Q = (lambdaZ * transpose(Ymat) * Ymat + rho * eye(N,N) + rho * onematrix);
    R = ((lambdaZ * transpose(Ymat) * (Ymat - E)) + rho * (onematrix + C) - onevector * transpose(lowDelta) - delta);
    
    %%% QA = R
    %%% A = inv(Q) * R; 
    %%% MATLAB TELLS ME NOT TO USE INV() BECAUSE IT IS SLOW AND MAY BE LESS ACCURATE
    %%% INSTEAD, IT WANTS ME TO USE THE \ OPERATION.
    
    %%% I TRUST YOU MATLAB
    A = Q\R;
    
    inJmatrix = A + delta / rho;
    
    %%% THIS WTHRESH FUNCTION IS A BUILT-IN SHRINKAGE-THRESHOLDING OPERATOR
    J = wthresh(inJmatrix,'s',1/rho);
    C = J - diag(diag(J));
    
    
    %%% lowDelta
    lowDelta = lowDelta + rho * (transpose(A) * onevector - onevector);
    
    %%% delta
    delta = delta + rho * (A - C);

%     
%     if norm(transpose(A) * onevector - onevector, inf) <= epsilon && norm(A-C, inf) <= epsilon && norm(A-prevA, inf) <= epsilon 
%         terminate = true;
%     end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           NORMALIZE C AND CONSTRUCT W          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I do not know if this is correct, 
%%% I am trying to normalize the columns of C 
%%% using ci = ci / norm(ci, inf)
%%% as stated in Algorithm 1 in the paper. (Page 5)

Cmat = zeros(N,N);
for i=1:N
    column = C(:,i);
    column = column / norm(column, inf);
    Cmat(:,i) = column;
end


%%% Construct the Wmat matrix
Wmat = abs(Cmat) + transpose(abs(Cmat));
%%% Create the figure
clim = [0 1];
imagesc(Wmat, clim)
