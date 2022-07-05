function outfun = spdaux_smootheval(input, h, x)

% auxiliary function that returns an evaluating functional given 'x'
% and returns a matrix that concerns a given input.

%% Preprocessing
p = input.size(1); % dimension
N = input.size(3); % number of data points
y = spdaux_equivariant(input); % (N-by-p^2)

%% Function Handle Creater
outfun = @spd_eval1d;

%% auxiliary function here
function value = spdaux_smooth_kernel(x,h)
value = exp(-(x^2)/(2*(h^2)))/h;
end


%% Then, What Function ?
    function output = spd_eval1d(xtgt)
        % 1. denominator, now a single number
        xdenom = 0;
        for j=1:N
            xdenom = xdenom + spdaux_smooth_kernel(xtgt-x(j),h);
        end
        % 2. compute the smoothed ones
        ytgt = zeros(1,p^2);
        for i=1:N
            ytgt = ytgt + y(i,:)*spdaux_smooth_kernel(x(i)-xtgt,h)/xdenom;
        end
        % 3. convert back
        y3d = reshape(ytgt, [p,p]);
        y3dsym = (y3d+(y3d'))/2;
        output = expm(y3dsym);
    end
end