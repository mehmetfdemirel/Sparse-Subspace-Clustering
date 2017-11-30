%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      MAKE DATA                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 360 * 640;
Drows = dim;
N1_size = 20;


Y1mat = reshape(rgb2gray(im2double(imread('00005.png'))), [dim,1]);
Y2mat = reshape(rgb2gray(im2double(imread('00010.png'))), [dim,1]);
Y3mat = reshape(rgb2gray(im2double(imread('00015.png'))), [dim,1]);
Y4mat = reshape(rgb2gray(im2double(imread('00020.png'))), [dim,1]);
Y5mat = reshape(rgb2gray(im2double(imread('00025.png'))), [dim,1]);
Y6mat = reshape(rgb2gray(im2double(imread('00030.png'))), [dim,1]);
Y7mat = reshape(rgb2gray(im2double(imread('00035.png'))), [dim,1]);
Y8mat = reshape(rgb2gray(im2double(imread('00040.png'))), [dim,1]);
Y9mat = reshape(rgb2gray(im2double(imread('00045.png'))), [dim,1]);
Y10mat = reshape(rgb2gray(im2double(imread('00050.png'))), [dim,1]);
Y11mat = reshape(rgb2gray(im2double(imread('00055.png'))), [dim,1]);
Y12mat = reshape(rgb2gray(im2double(imread('00060.png'))), [dim,1]);
Y13mat = reshape(rgb2gray(im2double(imread('00065.png'))), [dim,1]);
Y14mat = reshape(rgb2gray(im2double(imread('00070.png'))), [dim,1]);
Y15mat = reshape(rgb2gray(im2double(imread('00075.png'))), [dim,1]);
Y16mat = reshape(rgb2gray(im2double(imread('00080.png'))), [dim,1]);
Y17mat = reshape(rgb2gray(im2double(imread('00085.png'))), [dim,1]);
Y18mat = reshape(rgb2gray(im2double(imread('00090.png'))), [dim,1]);
Y19mat = reshape(rgb2gray(im2double(imread('00095.png'))), [dim,1]);
Y20mat = reshape(rgb2gray(im2double(imread('00100.png'))), [dim,1]);

Y21mat = reshape(rgb2gray(im2double(imread('01010.png'))), [dim,1]);
Y22mat = reshape(rgb2gray(im2double(imread('01020.png'))), [dim,1]);
Y23mat = reshape(rgb2gray(im2double(imread('01030.png'))), [dim,1]);
Y24mat = reshape(rgb2gray(im2double(imread('01040.png'))), [dim,1]);
Y25mat = reshape(rgb2gray(im2double(imread('01050.png'))), [dim,1]);
Y26mat = reshape(rgb2gray(im2double(imread('01060.png'))), [dim,1]);
Y27mat = reshape(rgb2gray(im2double(imread('01070.png'))), [dim,1]);
Y28mat = reshape(rgb2gray(im2double(imread('01080.png'))), [dim,1]);
Y29mat = reshape(rgb2gray(im2double(imread('01090.png'))), [dim,1]);
Y30mat = reshape(rgb2gray(im2double(imread('01100.png'))), [dim,1]);
Y31mat = reshape(rgb2gray(im2double(imread('01110.png'))), [dim,1]);
Y32mat = reshape(rgb2gray(im2double(imread('01120.png'))), [dim,1]);
Y33mat = reshape(rgb2gray(im2double(imread('01130.png'))), [dim,1]);
Y34mat = reshape(rgb2gray(im2double(imread('01140.png'))), [dim,1]);
Y35mat = reshape(rgb2gray(im2double(imread('01150.png'))), [dim,1]);
Y36mat = reshape(rgb2gray(im2double(imread('01160.png'))), [dim,1]);
Y37mat = reshape(rgb2gray(im2double(imread('01170.png'))), [dim,1]);
Y38mat = reshape(rgb2gray(im2double(imread('01180.png'))), [dim,1]);
Y39mat = reshape(rgb2gray(im2double(imread('01190.png'))), [dim,1]);
Y40mat = reshape(rgb2gray(im2double(imread('01200.png'))), [dim,1]);

Y41mat = reshape(rgb2gray(im2double(imread('02020.png'))), [dim,1]);
Y42mat = reshape(rgb2gray(im2double(imread('02040.png'))), [dim,1]);
Y43mat = reshape(rgb2gray(im2double(imread('02060.png'))), [dim,1]);
Y44mat = reshape(rgb2gray(im2double(imread('02080.png'))), [dim,1]);
Y45mat = reshape(rgb2gray(im2double(imread('02100.png'))), [dim,1]);
Y46mat = reshape(rgb2gray(im2double(imread('02120.png'))), [dim,1]);
Y47mat = reshape(rgb2gray(im2double(imread('02140.png'))), [dim,1]);
Y48mat = reshape(rgb2gray(im2double(imread('02160.png'))), [dim,1]);
Y49mat = reshape(rgb2gray(im2double(imread('02180.png'))), [dim,1]);
Y50mat = reshape(rgb2gray(im2double(imread('02200.png'))), [dim,1]);
Y51mat = reshape(rgb2gray(im2double(imread('02220.png'))), [dim,1]);
Y52mat = reshape(rgb2gray(im2double(imread('02240.png'))), [dim,1]);
Y53mat = reshape(rgb2gray(im2double(imread('02260.png'))), [dim,1]);
Y54mat = reshape(rgb2gray(im2double(imread('02280.png'))), [dim,1]);
Y55mat = reshape(rgb2gray(im2double(imread('02300.png'))), [dim,1]);
Y56mat = reshape(rgb2gray(im2double(imread('02320.png'))), [dim,1]);
Y57mat = reshape(rgb2gray(im2double(imread('02340.png'))), [dim,1]);
Y58mat = reshape(rgb2gray(im2double(imread('02360.png'))), [dim,1]);
Y59mat = reshape(rgb2gray(im2double(imread('02380.png'))), [dim,1]);
Y60mat = reshape(rgb2gray(im2double(imread('02400.png'))), [dim,1]);

Y61mat = reshape(rgb2gray(im2double(imread('03030.png'))), [dim,1]);
Y62mat = reshape(rgb2gray(im2double(imread('03060.png'))), [dim,1]);
Y63mat = reshape(rgb2gray(im2double(imread('03090.png'))), [dim,1]);
Y64mat = reshape(rgb2gray(im2double(imread('03120.png'))), [dim,1]);
Y65mat = reshape(rgb2gray(im2double(imread('03150.png'))), [dim,1]);
Y66mat = reshape(rgb2gray(im2double(imread('03180.png'))), [dim,1]);
Y67mat = reshape(rgb2gray(im2double(imread('03210.png'))), [dim,1]);
Y68mat = reshape(rgb2gray(im2double(imread('03240.png'))), [dim,1]);
Y69mat = reshape(rgb2gray(im2double(imread('03270.png'))), [dim,1]);
Y70mat = reshape(rgb2gray(im2double(imread('03300.png'))), [dim,1]);
Y71mat = reshape(rgb2gray(im2double(imread('03330.png'))), [dim,1]);
Y72mat = reshape(rgb2gray(im2double(imread('03360.png'))), [dim,1]);
Y73mat = reshape(rgb2gray(im2double(imread('03390.png'))), [dim,1]);
Y74mat = reshape(rgb2gray(im2double(imread('03420.png'))), [dim,1]);
Y75mat = reshape(rgb2gray(im2double(imread('03450.png'))), [dim,1]);
Y76mat = reshape(rgb2gray(im2double(imread('03480.png'))), [dim,1]);
Y77mat = reshape(rgb2gray(im2double(imread('03510.png'))), [dim,1]);
Y78mat = reshape(rgb2gray(im2double(imread('03540.png'))), [dim,1]);
Y79mat = reshape(rgb2gray(im2double(imread('03570.png'))), [dim,1]);
Y80mat = reshape(rgb2gray(im2double(imread('03600.png'))), [dim,1]);

Y81mat = reshape(rgb2gray(im2double(imread('04040.png'))), [dim,1]);
Y82mat = reshape(rgb2gray(im2double(imread('04080.png'))), [dim,1]);
Y83mat = reshape(rgb2gray(im2double(imread('04120.png'))), [dim,1]);
Y84mat = reshape(rgb2gray(im2double(imread('04160.png'))), [dim,1]);
Y85mat = reshape(rgb2gray(im2double(imread('04200.png'))), [dim,1]);
Y86mat = reshape(rgb2gray(im2double(imread('04240.png'))), [dim,1]);
Y87mat = reshape(rgb2gray(im2double(imread('04280.png'))), [dim,1]);
Y88mat = reshape(rgb2gray(im2double(imread('04320.png'))), [dim,1]);
Y89mat = reshape(rgb2gray(im2double(imread('04360.png'))), [dim,1]);
Y90mat = reshape(rgb2gray(im2double(imread('04400.png'))), [dim,1]);
Y91mat = reshape(rgb2gray(im2double(imread('04440.png'))), [dim,1]);
Y92mat = reshape(rgb2gray(im2double(imread('04480.png'))), [dim,1]);
Y93mat = reshape(rgb2gray(im2double(imread('04520.png'))), [dim,1]);
Y94mat = reshape(rgb2gray(im2double(imread('04560.png'))), [dim,1]);
Y95mat = reshape(rgb2gray(im2double(imread('04600.png'))), [dim,1]);
Y96mat = reshape(rgb2gray(im2double(imread('04640.png'))), [dim,1]);
Y97mat = reshape(rgb2gray(im2double(imread('04680.png'))), [dim,1]);
Y98mat = reshape(rgb2gray(im2double(imread('04720.png'))), [dim,1]);
Y99mat = reshape(rgb2gray(im2double(imread('04760.png'))), [dim,1]);
Y100mat = reshape(rgb2gray(im2double(imread('04800.png'))), [dim,1]);

YY1mat = [Y1mat Y2mat Y3mat Y4mat Y5mat Y6mat Y7mat Y8mat Y9mat Y10mat Y11mat Y12mat Y13mat Y14mat Y15mat Y16mat Y17mat Y18mat Y19mat Y20mat];
YY2mat = [Y21mat Y22mat Y23mat Y24mat Y25mat Y26mat Y27mat Y28mat Y29mat Y30mat Y31mat Y32mat Y33mat Y34mat Y35mat Y36mat Y37mat Y38mat Y39mat Y40mat];
YY3mat = [Y41mat Y42mat Y43mat Y44mat Y45mat Y46mat Y47mat Y48mat Y49mat Y50mat Y51mat Y52mat Y53mat Y54mat Y55mat Y56mat Y57mat Y58mat Y59mat Y60mat];
YY4mat = [Y61mat Y62mat Y63mat Y64mat Y65mat Y66mat Y67mat Y68mat Y69mat Y70mat Y71mat Y72mat Y73mat Y74mat Y75mat Y76mat Y77mat Y78mat Y79mat Y80mat];
YY5mat = [Y81mat Y82mat Y83mat Y84mat Y85mat Y86mat Y87mat Y88mat Y89mat Y90mat Y91mat Y92mat Y93mat Y94mat Y95mat Y96mat Y97mat Y98mat Y99mat Y100mat];
Ymat = [YY1mat YY2mat YY3mat YY4mat YY5mat];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    INITILIZATION                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = Drows;
N = N1_size * 5;
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%