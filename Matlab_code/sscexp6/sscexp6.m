%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      MAKE DATA                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 360 * 640;
Drows = dim;
N1_size = 7;


Y1mat = reshape(rgb2gray(im2double(imread('00001.png'))), [dim,1]);
Y2mat = reshape(rgb2gray(im2double(imread('00002.png'))), [dim,1]);
Y3mat = reshape(rgb2gray(im2double(imread('00003.png'))), [dim,1]);
Y4mat = reshape(rgb2gray(im2double(imread('00005.png'))), [dim,1]);
Y5mat = reshape(rgb2gray(im2double(imread('00008.png'))), [dim,1]);
Y6mat = reshape(rgb2gray(im2double(imread('00013.png'))), [dim,1]);
Y7mat = reshape(rgb2gray(im2double(imread('00021.png'))), [dim,1]);

Y8mat = reshape(rgb2gray(im2double(imread('00034.png'))), [dim,1]);
Y9mat = reshape(rgb2gray(im2double(imread('00055.png'))), [dim,1]);
Y10mat = reshape(rgb2gray(im2double(imread('00089.png'))), [dim,1]);
Y11mat = reshape(rgb2gray(im2double(imread('00144.png'))), [dim,1]);
Y12mat = reshape(rgb2gray(im2double(imread('00233.png'))), [dim,1]);
Y13mat = reshape(rgb2gray(im2double(imread('00377.png'))), [dim,1]);
Y14mat = reshape(rgb2gray(im2double(imread('00610.png'))), [dim,1]);

Y15mat = reshape(rgb2gray(im2double(imread('00987.png'))), [dim,1]);
Y16mat = reshape(rgb2gray(im2double(imread('01597.png'))), [dim,1]);
Y17mat = reshape(rgb2gray(im2double(imread('02584.png'))), [dim,1]);
Y18mat = reshape(rgb2gray(im2double(imread('04181.png'))), [dim,1]);
Y19mat = reshape(rgb2gray(im2double(imread('06765.png'))), [dim,1]);
Y20mat = reshape(rgb2gray(im2double(imread('10946.png'))), [dim,1]);
Y21mat = reshape(rgb2gray(im2double(imread('15691.png'))), [dim,1]);


YY1mat = [Y1mat Y2mat Y3mat Y4mat Y5mat Y6mat Y7mat];
YY2mat = [Y8mat Y9mat Y10mat Y11mat Y12mat Y13mat Y14mat];
YY3mat= [Y15mat Y16mat Y17mat Y18mat Y19mat Y20mat Y21mat];

Ymat = [YY1mat YY2mat YY3mat];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    INITILIZATION                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = Drows;
N = N1_size * 3;
maxIter = 200;

C = zeros(N,N);
A = zeros(N,N);

E = zeros(D,N);

lowDelta = zeros(N,1);
delta = zeros(N,N);
rho = 100;
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
    %%% prevA = A;
    Q = (lambdaZ * transpose(Ymat) * Ymat + rho * eye(N,N) + rho * onematrix);
    R = ((lambdaZ * transpose(Ymat) * (Ymat - E)) + rho * (onematrix + C) - onevector * transpose(lowDelta) - delta);
    
    %%% QA = R
    %%%A = inv(Q) * R; 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               NORMALIZE C AND CONSTRUCT W               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

sum(C(:,3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%