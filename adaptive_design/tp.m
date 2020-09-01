function [y] = tp(x)
l = size(x, 1);
mean = zeros(l, 1);
y = random('tlocationscale', mean, 1, 3, l, 1);
end