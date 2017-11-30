%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      MAKE DATA                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 360 * 640;
Drows = dim;
N1_size = 4;


Y1mat = reshape(rgb2gray(im2double(imread('00900.png'))), [dim,1]);
Y2mat = reshape(rgb2gray(im2double(imread('00925.png'))), [dim,1]);
Y3mat = reshape(rgb2gray(im2double(imread('00950.png'))), [dim,1]);
Y4mat = reshape(rgb2gray(im2double(imread('00975.png'))), [dim,1]);
Y5mat = reshape(rgb2gray(im2double(imread('12000.png'))), [dim,1]);
Y6mat = reshape(rgb2gray(im2double(imread('12025.png'))), [dim,1]);
Y7mat = reshape(rgb2gray(im2double(imread('12050.png'))), [dim,1]);
Y8mat = reshape(rgb2gray(im2double(imread('12075.png'))), [dim,1]);


%YY1mat = [Y1mat Y2mat Y3mat Y4mat Y5mat Y6mat Y7mat Y8mat Y9mat Y10mat Y11mat Y12mat Y13mat Y14mat Y15mat Y16mat Y17mat Y18mat Y19mat Y20mat Y21mat Y22mat Y23mat Y24mat Y25mat Y26mat Y27mat Y28mat Y29mat Y30mat Y31mat Y32mat Y33mat Y34mat Y35mat Y36mat Y37mat Y38mat Y39mat Y40mat Y41mat Y42mat Y43mat Y44mat Y45mat Y46mat Y47mat Y48mat Y49mat Y50mat Y51mat Y52mat Y53mat Y54mat Y55mat Y56mat Y57mat Y58mat Y59mat Y60mat Y61mat Y62mat Y63mat Y64mat Y65mat Y66mat Y67mat Y68mat Y69mat Y70mat Y71mat Y72mat Y73mat Y74mat Y75mat Y76mat Y77mat Y78mat Y79mat Y80mat Y81mat Y82mat Y83mat Y84mat Y85mat Y86mat Y87mat Y88mat Y89mat Y90mat Y91mat Y92mat Y93mat Y94mat Y95mat Y96mat Y97mat Y98mat Y99mat Y100mat Y101mat Y102mat Y103mat Y104mat Y105mat Y106mat Y107mat Y108mat Y109mat Y110mat Y111mat Y112mat Y113mat Y114mat Y115mat Y116mat Y117mat Y118mat Y119mat Y120mat Y121mat Y122mat Y123mat Y124mat Y125mat Y126mat Y127mat Y128mat Y129mat Y130mat Y131mat Y132mat Y133mat Y134mat Y135mat Y136mat Y137mat Y138mat Y139mat Y140mat Y141mat Y142mat Y143mat Y144mat Y145mat Y146mat Y147mat Y148mat Y149mat Y150mat];
YY1mat = [Y1mat Y2mat Y3mat Y4mat];
YY2mat = [Y5mat Y6mat Y7mat Y8mat];

Ymat = [YY1mat YY2mat];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    INITILIZATION                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = Drows;
N = N1_size * 2;
maxIter = 200;

C = zeros(N,N);
A = zeros(N,N);

E = zeros(D,N);

lowDelta = zeros(N,1);
delta = zeros(N,N);
rho = 1000;
lambdaZ = 1;
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


sum(C(:,3))