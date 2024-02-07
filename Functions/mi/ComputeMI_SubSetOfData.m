function Vargout=ComputeMI_SubSetOfData(xs,ys)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE MI FROM SCRATCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
The goal is to code  p^(x)=1N∑i=1Nδ(x−x(i),h), eq 17 in Peng 2005,
which is the joint probability of x and y. 
This required the following steps:
    1. Define δ 
        1.1 Define z, h, sigma
        1.2 Check for peer review - is sigma positive and (semi)definite?
        1.3 Compute δ 
    2. Compute the approximate density function
    3. Compare to the MI results using Peng's functions
%}


points=[xs ys];


for ii=1:length(xs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Compute δ
% STEP 1.1: Define z, h, and Sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
The equation governing δ is 
δ(z,h)=exp(−zTΣ−1z2h2) / {(2π)d/2hd|Σ|1/2}
As inputs, this needs z, h, and Sigma
%}

% Getting z, which is x-xi - extract the i-th data point and convert it to a column vector
xi = xs(:,ii); 
z=xs-xi;

% Define sigma, which is the covariance of z
Sigma=cov(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1.2
% CRUCIAL FOR PEER REVIEW - CHECK TO SEE IF THE COVARIANCE MATRIX IS 
% POSITIVE AND (SEMI)DEFINITE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
The reviewer is concerned it is not because the number of time points
(42) is comparable to the number of time-series (100), 
the covariance matrix will not be positive semi-definite as expected, 
and a regularisation procedure has to be implemented. 

This website - 
https://math.stackexchange.com/questions/2533610/how-do-you-recognize-a-positive-semidefinite-matrix
Gave three examples to compute whether a matrix is positive and definite.
One that I understand is that if the matrix is (1) symmetric and 
(2) and all eigenvalues are positive, then the matrix is positive definite.

And this thread gave me the code to do that - 
https://uk.mathworks.com/matlabcentral/answers/589435-how-do-i-create-a-function-script-to-check-the-positive-definiteness-of-a-a-square-matrix-of-any-siz
pos_def = isequal(Sigma,Sigma') && all( eig(Sigma)>0 );
%}

pos_def = isequal(Sigma,flip(Sigma)) && all( eig(Sigma)>0 );

% Store when the cov matrix is and isn't positive & semi definite
% Column 1 will be all the "good" times, col 2 will be "bad"
if isequal(pos_def,1)
    GoodAndBad(ii,1:2)=[1,0];
else
    GoodAndBad(ii,1:2)=[0,1];
end

%{
Define h, which is window width
If we want to compute h that scales with the variance around the point,
I've taken code from the below repo, and it's helped
https://gist.github.com/zahlenteufel/0ea490a1b28b7cc62492
To justify what h is to the reviewer, I would say that their previous
work (and cite the paper I used for ChatGPT) defined the window around h
as a Gaussian
%}

M = ii;
for i = 1:M
for j= 1:M
    pXY(i,j) = 0;
    
    Hmax = M;
    Hmin = 0.001;
    
    while Hmax - Hmin > 0.1
        H = (Hmax+Hmin) / 2;
        inside = sum(pointsdist(points, [xs(i) ys(j)]') < H);
        if inside > sqrt(M)
            Hmax = H;
        else 
            Hmin = H;
        end
    end
    
    h = H;

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1.3: Compute δ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
Calculates the value of Parzen window function, 
equation 18 in Peng 2005 paper
https://ieeexplore-ieee-org.surrey.idm.oclc.org/document/1453511/
δ(z,h)=exp(−zTΣ−1z2h2) / {(2π)d/2hd|Σ|1/2}
%}

d = length(z); % dimensionality of the problem

% calculate the denominator
denom = (2*pi)^(d/2) * h^d * sqrt(det(Sigma));

% calculate the numerator
num = exp(-z' * inv(Sigma) * z / (2 * h^2));

% calculate the result 
DELTA{ii,1} = num / denom;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Compute the approximate density function
% by adding δ, or delta, back to desity estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Given N samples of a variable x, the approximate density function p^(x)
computed via Parzen window. It has the following form, which is eq 17:
p^(x)=1N∑i=1Nδ(x−x(i),h) 
where δ(.) is the Parzen window function as above
x(i) is the ith sample
h is the window width.

The inline function I was recommended to compute a multivariate Gaussian 
kernel function with known bandwidth doesn't work all in one, so I've
broken it down
kernel = @(z, h) (2*pi*h^2)^(-length(z)/2) * exp(-0.5*z'*inv(h^2*eye(length(z)))*z); %
%}

result = 0;
N = size(DELTA{ii,1}, 1);

for i = 1:N

    % result = result + kernel(ii - DELTA{ii,1}(i,:), h);
    firstTerm=(2*pi*h^2)^(-length(z)/2);
    secondTermPart1=-0.5*z';
    secondTermPart2=inv(h^2*eye(length(z)));
    secondTermPart3=secondTermPart2.*z;
    secondTerm=exp(secondTermPart1.*secondTermPart3);
    kernel=firstTerm*secondTerm;
    result=result+kernel;

end

Result = result / N; % divide by N to get the average
tmpVargout = Result(1,1); % divide by N to get the average

end
Vargout=mean(tmpVargout)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Compare scratch miFC with that computed by Peng's algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute miFC using Peng's function
Mutinf=mutualinfo(xs,ys);
end