%%%%%%%
% CODE BELOW WAS TAKEN FROM
% https://gist.github.com/zahlenteufel/0ea490a1b28b7cc62492


% Generate the points
clf
subplot(2,1,1);
axis equal
mu1 = [5 4]';
sigma1 = [2 -0.8; -0.8 1];
mu2 = [0 1]';
sigma2 = [1 -0.06;-0.06 0.5];
chol(sigma1);
chol(sigma2); % verifica definida positiva
P1 = 0.3;
P2 = 1 - P1;
n = 500;
n1 = sum(rand(1, n) <= P1);
n2 = n - n1;
x1 = bsxfun(@plus, sigma1 * randn(2, n1), mu1);
x2 = bsxfun(@plus, sigma2 * randn(2, n2), mu2);
points = [x1 x2];

plot(x1(1,:), x1(2,:), 'r+');
hold on
plot(x2(1,:), x2(2,:), 'b+');
xlim([0 7]);
ylim([0 7]);

M = 50;
xs = linspace(0, 7, M);
ys = linspace(0, 7, M);

pXY = zeros(M);
[X Y] = meshgrid(xs,ys);

%
%h = 0.25;
hold on
for i = 1:M
for j= 1:M
    pXY(i,j) = 0;
    
    Hmax = 50;
    Hmin = 0.001;
    
    while Hmax - Hmin > 0.1
        H = (Hmax+Hmin) / 2;
        inside = sum(pointsdist(points, [xs(i) ys(j)]') < H);
        if inside > sqrt(n)
            Hmax = H;
        else 
            Hmin = H;
        end
    end
    
    h = H;
        
    %circle(xs(i), ys(j), h);
    
    for k=1:n             
        dist = norm([xs(i) ys(j)]' - points(:,k));
        pXY(i,j) = pXY(i,j) + normpdf(dist/h,0,1);
    end
end
end

subplot(2,1,2);
surf(X, Y, -pXY');
title(['h = ', num2str(h)]);
xlim([0 7]);
ylim([0 7]);
axis square

fprintf('n = %d, sqrt(n) = %f \n', n, sqrt(n));
subplot(2,1,1);
hold off

plot(x1(1,:), x1(2,:), 'r+');
hold on
plot(x2(1,:), x2(2,:), 'b+');
xlim([0 7]);
ylim([0 7]);
axis square
box off
% 
while 1
    [mx,my] = ginput(1);
     for H = 0.01:0.2:10
            inside = sum(pointsdist(points, [mx my]') < H);
            if inside > sqrt(n)
                h = H;
                break;
            end
     end
       circle(mx, my, h);
end


% 
function d = pointsdist(points, x)
% calculates the distance from each point to x
difs = bsxfun(@minus, points, x);
d = sqrt(sum(difs.^2,1));
end