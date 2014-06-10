function [d,gamma] = discrimVectors(x, K, c)
%start

N = [size(x{1}, 1) size(x{2}, 1)];

featureCount = size(x{1}(1, :), 2);

%estimated mean of class
mi = [
    sum(x{1}, 1)/N(1)
    sum(x{2}, 1)/N(2)
];

%diffeerence in the estimated means
delta = transpose(mi(1, :) - mi(2, :));

%within-class scatter for class
W = {0;0};

for j=1:N(1)    
    e = (x{1}(j, :) - mi(1, :)) * transpose(x{1}(j, :) - mi(1, :));
    
    W{1} = W{1} + e;
end

for j=1:N(2)
    e = (x{2}(j, :) - mi(2, :))* transpose(x{2}(j, :) - mi(2, :));
    
    W{2} = W{2} + e;
end

Wsum = W{1} + W{2};

A = double(c * W{1} + (1 - c) * W{2});

%compute d1 and S1
d = zeros(K + 1, featureCount);
S = [];

alfa1 = sqrt(inv(transpose(delta) * (1/A)^2 * delta));

d(1, 1:featureCount) = normalize(transpose(inv(A) * delta));

s11 = transpose(d(:, 1)) * inv(A) * d(:, 1); %w art nie ma tej transpozycji

n = 1;
S_inv = {};
S_inv{1} = 1/double(s11);

s = zeros(K, K);
s(1, 1) = s11;

gamma = zeros(K, 1);

%loop
while(1==1)    
    n = n + 1;

    %compute d_n           
    
    z = eye(n - 1); 
    z = z(:, 1) / alfa1;           
    
    tmp = transpose((d(1:n-1, 1:featureCount))) ...
        * S_inv{n - 1}...
        * z;
   
    new_d = double(inv(A)) * (delta - tmp);
    
    d(n, 1:featureCount) = normalize(new_d);
    
    %compute gamma_n                          
    gamma(n) = ((d(n, 1:featureCount) * delta)^2) ...
    /(d(n, 1:featureCount) * double(A) * transpose(d(n, 1:featureCount)));
    
    if n == K + 1
        break 
    end

    %compute S_n
    for i=1:n
        s(i, n) = (d(i, :)) * inv(A) * transpose(d(n, :));
    end
    
    for i=1:n-1
        s(n, i) = (d(n, :)) * inv(A) * transpose(d(i, :));
    end       
    
    y = transpose(s(n, 1:n-1));
    
    c_n = s(n, n) - transpose(y) * S_inv{n - 1} * y;                
    
    tmp_a = c_n * S_inv{n - 1} + S_inv{n - 1} * y * transpose(y) * S_inv{n - 1};
    tmp_b = -S_inv{n - 1} * y;
    tmp_c = -transpose(y) * S_inv{n - 1};
    tmp_d = 1;
    
    S_inv{n} = [tmp_a tmp_b; tmp_c tmp_d]/c_n;
end

