%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      MAKE DATA                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 360 * 640;
Drows = dim;
N1_size = 30;


Y1mat = reshape(rgb2gray(im2double(imread('06000.png'))), [dim,1]);
Y2mat = reshape(rgb2gray(im2double(imread('06024.png'))), [dim,1]);
Y3mat = reshape(rgb2gray(im2double(imread('06048.png'))), [dim,1]);
Y4mat = reshape(rgb2gray(im2double(imread('06072.png'))), [dim,1]);
Y5mat = reshape(rgb2gray(im2double(imread('06096.png'))), [dim,1]);
Y6mat = reshape(rgb2gray(im2double(imread('06120.png'))), [dim,1]);
Y7mat = reshape(rgb2gray(im2double(imread('06144.png'))), [dim,1]);
Y8mat = reshape(rgb2gray(im2double(imread('06168.png'))), [dim,1]);
Y9mat = reshape(rgb2gray(im2double(imread('06192.png'))), [dim,1]);
Y10mat = reshape(rgb2gray(im2double(imread('06216.png'))), [dim,1]);
Y11mat = reshape(rgb2gray(im2double(imread('06240.png'))), [dim,1]);
Y12mat = reshape(rgb2gray(im2double(imread('06264.png'))), [dim,1]);
Y13mat = reshape(rgb2gray(im2double(imread('06288.png'))), [dim,1]);
Y14mat = reshape(rgb2gray(im2double(imread('06312.png'))), [dim,1]);
Y15mat = reshape(rgb2gray(im2double(imread('06336.png'))), [dim,1]);
Y16mat = reshape(rgb2gray(im2double(imread('06360.png'))), [dim,1]);
Y17mat = reshape(rgb2gray(im2double(imread('06384.png'))), [dim,1]);
Y18mat = reshape(rgb2gray(im2double(imread('06408.png'))), [dim,1]);
Y19mat = reshape(rgb2gray(im2double(imread('06432.png'))), [dim,1]);
Y20mat = reshape(rgb2gray(im2double(imread('06456.png'))), [dim,1]);
Y21mat = reshape(rgb2gray(im2double(imread('06480.png'))), [dim,1]);
Y22mat = reshape(rgb2gray(im2double(imread('06504.png'))), [dim,1]);
Y23mat = reshape(rgb2gray(im2double(imread('06528.png'))), [dim,1]);
Y24mat = reshape(rgb2gray(im2double(imread('06552.png'))), [dim,1]);
Y25mat = reshape(rgb2gray(im2double(imread('06576.png'))), [dim,1]);
Y26mat = reshape(rgb2gray(im2double(imread('06600.png'))), [dim,1]);
Y27mat = reshape(rgb2gray(im2double(imread('06624.png'))), [dim,1]);
Y28mat = reshape(rgb2gray(im2double(imread('06648.png'))), [dim,1]);
Y29mat = reshape(rgb2gray(im2double(imread('06672.png'))), [dim,1]);
Y30mat = reshape(rgb2gray(im2double(imread('06696.png'))), [dim,1]);


Y31mat = reshape(rgb2gray(im2double(imread('08000.png'))), [dim,1]);
Y32mat = reshape(rgb2gray(im2double(imread('08024.png'))), [dim,1]);
Y33mat = reshape(rgb2gray(im2double(imread('08048.png'))), [dim,1]);
Y34mat = reshape(rgb2gray(im2double(imread('08072.png'))), [dim,1]);
Y35mat = reshape(rgb2gray(im2double(imread('08096.png'))), [dim,1]);
Y36mat = reshape(rgb2gray(im2double(imread('08120.png'))), [dim,1]);
Y37mat = reshape(rgb2gray(im2double(imread('08144.png'))), [dim,1]);
Y38mat = reshape(rgb2gray(im2double(imread('08168.png'))), [dim,1]);
Y39mat = reshape(rgb2gray(im2double(imread('08192.png'))), [dim,1]);
Y40mat = reshape(rgb2gray(im2double(imread('08216.png'))), [dim,1]);
Y41mat = reshape(rgb2gray(im2double(imread('08240.png'))), [dim,1]);
Y42mat = reshape(rgb2gray(im2double(imread('08264.png'))), [dim,1]);
Y43mat = reshape(rgb2gray(im2double(imread('08288.png'))), [dim,1]);
Y44mat = reshape(rgb2gray(im2double(imread('08312.png'))), [dim,1]);
Y45mat = reshape(rgb2gray(im2double(imread('08336.png'))), [dim,1]);
Y46mat = reshape(rgb2gray(im2double(imread('08360.png'))), [dim,1]);
Y47mat = reshape(rgb2gray(im2double(imread('08384.png'))), [dim,1]);
Y48mat = reshape(rgb2gray(im2double(imread('08408.png'))), [dim,1]);
Y49mat = reshape(rgb2gray(im2double(imread('08432.png'))), [dim,1]);
Y50mat = reshape(rgb2gray(im2double(imread('08456.png'))), [dim,1]);
Y51mat = reshape(rgb2gray(im2double(imread('08480.png'))), [dim,1]);
Y52mat = reshape(rgb2gray(im2double(imread('08504.png'))), [dim,1]);
Y53mat = reshape(rgb2gray(im2double(imread('08528.png'))), [dim,1]);
Y54mat = reshape(rgb2gray(im2double(imread('08552.png'))), [dim,1]);
Y55mat = reshape(rgb2gray(im2double(imread('08576.png'))), [dim,1]);
Y56mat = reshape(rgb2gray(im2double(imread('08600.png'))), [dim,1]);
Y57mat = reshape(rgb2gray(im2double(imread('08624.png'))), [dim,1]);
Y58mat = reshape(rgb2gray(im2double(imread('08648.png'))), [dim,1]);
Y59mat = reshape(rgb2gray(im2double(imread('08672.png'))), [dim,1]);
Y60mat = reshape(rgb2gray(im2double(imread('08696.png'))), [dim,1]);


Y61mat = reshape(rgb2gray(im2double(imread('10000.png'))), [dim,1]);
Y62mat = reshape(rgb2gray(im2double(imread('10024.png'))), [dim,1]);
Y63mat = reshape(rgb2gray(im2double(imread('10048.png'))), [dim,1]);
Y64mat = reshape(rgb2gray(im2double(imread('10072.png'))), [dim,1]);
Y65mat = reshape(rgb2gray(im2double(imread('10096.png'))), [dim,1]);
Y66mat = reshape(rgb2gray(im2double(imread('10120.png'))), [dim,1]);
Y67mat = reshape(rgb2gray(im2double(imread('10144.png'))), [dim,1]);
Y68mat = reshape(rgb2gray(im2double(imread('10168.png'))), [dim,1]);
Y69mat = reshape(rgb2gray(im2double(imread('10192.png'))), [dim,1]);
Y70mat = reshape(rgb2gray(im2double(imread('10216.png'))), [dim,1]);
Y71mat = reshape(rgb2gray(im2double(imread('10240.png'))), [dim,1]);
Y72mat = reshape(rgb2gray(im2double(imread('10264.png'))), [dim,1]);
Y73mat = reshape(rgb2gray(im2double(imread('10288.png'))), [dim,1]);
Y74mat = reshape(rgb2gray(im2double(imread('10312.png'))), [dim,1]);
Y75mat = reshape(rgb2gray(im2double(imread('10336.png'))), [dim,1]);
Y76mat = reshape(rgb2gray(im2double(imread('10360.png'))), [dim,1]);
Y77mat = reshape(rgb2gray(im2double(imread('10384.png'))), [dim,1]);
Y78mat = reshape(rgb2gray(im2double(imread('10408.png'))), [dim,1]);
Y79mat = reshape(rgb2gray(im2double(imread('10432.png'))), [dim,1]);
Y80mat = reshape(rgb2gray(im2double(imread('10456.png'))), [dim,1]);
Y81mat = reshape(rgb2gray(im2double(imread('10480.png'))), [dim,1]);
Y82mat = reshape(rgb2gray(im2double(imread('10504.png'))), [dim,1]);
Y83mat = reshape(rgb2gray(im2double(imread('10528.png'))), [dim,1]);
Y84mat = reshape(rgb2gray(im2double(imread('10552.png'))), [dim,1]);
Y85mat = reshape(rgb2gray(im2double(imread('10576.png'))), [dim,1]);
Y86mat = reshape(rgb2gray(im2double(imread('10600.png'))), [dim,1]);
Y87mat = reshape(rgb2gray(im2double(imread('10624.png'))), [dim,1]);
Y88mat = reshape(rgb2gray(im2double(imread('10648.png'))), [dim,1]);
Y89mat = reshape(rgb2gray(im2double(imread('10672.png'))), [dim,1]);
Y90mat = reshape(rgb2gray(im2double(imread('10696.png'))), [dim,1]);


YY1mat = [Y1mat Y2mat Y3mat Y4mat Y5mat Y6mat Y7mat Y8mat Y9mat Y10mat Y11mat Y12mat Y13mat Y14mat Y15mat Y16mat Y17mat Y18mat Y19mat Y20mat Y21mat Y22mat Y23mat Y24mat Y25mat Y26mat Y27mat Y28mat Y29mat Y30mat];
YY2mat = [Y31mat Y32mat Y33mat Y34mat Y35mat Y36mat Y37mat Y38mat Y39mat Y40mat Y41mat Y42mat Y43mat Y44mat Y45mat Y46mat Y47mat Y48mat Y49mat Y50mat Y51mat Y52mat Y53mat Y54mat Y55mat Y56mat Y57mat Y58mat Y59mat Y60mat];
YY3mat= [Y61mat Y62mat Y63mat Y64mat Y65mat Y66mat Y67mat Y68mat Y69mat Y70mat Y71mat Y72mat Y73mat Y74mat Y75mat Y76mat Y77mat Y78mat Y79mat Y80mat Y81mat Y82mat Y83mat Y84mat Y85mat Y86mat Y87mat Y88mat Y89mat Y90mat];

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

sum(C(:,4))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%