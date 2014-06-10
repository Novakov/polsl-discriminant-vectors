x = {};
x{1} = [2 5;2 6;3 5;3 6]; %class 1
%x{1} = [6 7;6 8;7 7;7 8]; %class 2
x{2} = [6 2;6 3;7 2;7 3]; %class 2
%x{1} = [-2 2;-5 0;-2 3; -1 4];
%x{2} = [3 3;5 1;3 -2;1 -5];


c = 0.5;

K = 6; %discriminant vector count

%start

N = [size(x{1}, 1) size(x{2}, 1)];

%estimated mean of class
mi = [
    sum(x{1}, 1)/N(1)
    sum(x{2}, 1)/N(2)
];

%diffeerence in the estimated means
delta = ([
        mi(1, 1) - mi(2, 1);
        mi(1, 2) - mi(2, 2);
]);

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
d = zeros(K + 1, 2);
S = [];

alfa1 = sqrt(inv(transpose(delta) * (1/A)^2 * delta));

d(1, 1:2) = transpose(alfa1 * inv(A) * delta); %2 do jakiejœ zmiennej

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
    
    z = eye(n - 1); %TODO do zmiennej
    z = z(:, 1) / alfa1;           
    
    tmp = transpose((d(1:n-1, 1:2))) ...
        * S_inv{n - 1}...
        * z;
   
    new_d = double(inv(A)) * (delta - tmp);
    
    alfa = sqrt(new_d(1)^2 + new_d(2)^2);
    
    d(n, 1:2) = new_d/alfa;
    
    %compute gamma_n                          
    gamma(n) = ((d(n, 1:2) * delta)^2) ...
    /(d(n, 1:2) * double(A) * transpose(d(n, 1:2)));
    
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

d
gamma

bar(abs(transpose(d)))
