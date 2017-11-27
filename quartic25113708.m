function rts = quartic25113708(C)

% |real| <= |a| + |b| + |c| + 1 -- covered in discussion
% There also exists two real, one positive and one negative
max_root_value = abs(C(1)) + abs(C(2)) + abs(C(3)) + 1;
min_root_value = -1 * max_root_value;

% Since two real are guaranteed, I will find the two real
% and then use the quadratic formula to find the other real
% or imaginary.

% Newton's Method
rts = zeros(1, 4);
max_iterations = 100;
coefficients = [1, C(1), C(2), C(3), -1];
deriv_coefficients = [0, 4, 3 * C(1), 2 * C(2), C(3)];

% Finding root starting from min_root_value
p_zero = min_root_value;

for i = 1:max_iterations
    f = coefficients(1) * p_zero.^4 + coefficients(2) * p_zero.^3 + coefficients(3) * p_zero.^2 + coefficients(4) * p_zero + coefficients(5);
    f_prime = deriv_coefficients(1) * p_zero.^4 + deriv_coefficients(2) * p_zero.^3 + deriv_coefficients(3) * p_zero.^2 + deriv_coefficients(4) * p_zero + deriv_coefficients(5);
    
    p = p_zero - (f ./ f_prime);
    
    % Check tolerance
    if abs(p - p_zero) < 10.^-5
        rts(1) = p;
        break
    else
        p_zero = p;
    end
end
rts(1) = p;

% Finding root starting from max_root_value
p_zero = max_root_value;

for i = 1:max_iterations
    f = coefficients(1) * p_zero.^4 + coefficients(2) * p_zero.^3 + coefficients(3) * p_zero.^2 + coefficients(4) * p_zero + coefficients(5);
    f_prime = deriv_coefficients(1) * p_zero.^4 + deriv_coefficients(2) * p_zero.^3 + deriv_coefficients(3) * p_zero.^2 + deriv_coefficients(4) * p_zero + deriv_coefficients(5);
    
    p = p_zero - (f ./ f_prime);
    
    % Check tolerance
    if abs(p - p_zero) < 10.^-5
        rts(2) = p;
        break
    else
        p_zero = p;
    end
end
rts(2) = p;

% Factor polynomial based on first two real - synthetic division?
[q, r] = deconv(deconv(coefficients, [1, -1 * rts(1)]), [1, -1 * rts(2)]);

real = (-1 * q(2)) / (2 * q(1));
potential_imaginary = sqrt(abs((q(2).^2) - (4 * q(1) * q(3)))) / (2 * q(1));
% imaginary
if (q(2).^2) - (4 * q(1) * q(3)) < 0
    rts_3 = complex(real,  potential_imaginary);
    rts_4 = complex(real, -1 * potential_imaginary);
else
    rts_3 = real + potential_imaginary;
    rts_4 = real - potential_imaginary;
end

% Refining guess for rts 3 with Newton's
p_zero = rts_3;

for i = 1:max_iterations
    f = coefficients(1) * p_zero.^4 + coefficients(2) * p_zero.^3 + coefficients(3) * p_zero.^2 + coefficients(4) * p_zero + coefficients(5);
    f_prime = deriv_coefficients(1) * p_zero.^4 + deriv_coefficients(2) * p_zero.^3 + deriv_coefficients(3) * p_zero.^2 + deriv_coefficients(4) * p_zero + deriv_coefficients(5);
    
    p = p_zero - (f ./ f_prime);
    
    % Check tolerance
    if abs(p - p_zero) < 10.^-5
        rts(3) = p;
        break
    else
        p_zero = p;
    end
end
rts(3) = p;

% Refining guess for rts 4 with Newton's
p_zero = rts_4;

for i = 1:max_iterations
    f = coefficients(1) * p_zero.^4 + coefficients(2) * p_zero.^3 + coefficients(3) * p_zero.^2 + coefficients(4) * p_zero + coefficients(5);
    f_prime = deriv_coefficients(1) * p_zero.^4 + deriv_coefficients(2) * p_zero.^3 + deriv_coefficients(3) * p_zero.^2 + deriv_coefficients(4) * p_zero + deriv_coefficients(5);
    
    p = p_zero - (f ./ f_prime);
    
    % Check tolerance
    if abs(p - p_zero) < 10.^-5
        rts(4) = p;
        break
    else
        p_zero = p;
    end
end
rts(4) = p;

end



