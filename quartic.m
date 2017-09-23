function [root1, root2, root3, root4] = quartic(C)

% |real roots| <= |a| + |b| + |c| + 1 -- covered in discussion
% There also exists two real roots, one positive and one negative
max_root_value = abs(C(1)) + abs(C(2)) + abs(C(3)) + 1
min_root_value = -1 * max_root_value

% Since two real roots are guaranteed, I will find the two real roots
% and then use the quadratic formula to find the other real roots
% or imaginary roots.

% Newton's Method

coefficients = [1, C(1), C(2), C(3), -1]
roots = [];
starting_value = min_root_value

% Finding two real roots
while length(roots) < 2
    derivative_coefficients = [0, coefficients(1) * 4, coefficients(2) * 3, coefficients(3) * 2, coefficients(4) * 1]
    x = starting_value
    x_n1 = x ...
           - ((coefficients(1) * x^4 + coefficients(2) * x^3 + coefficients(3) * x^2 + coefficients(4) * x + coefficients(5)) ...
              \(derivative_coefficients(
    
    

end

