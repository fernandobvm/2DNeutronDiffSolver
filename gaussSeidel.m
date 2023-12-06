function [x,tol_reached] = gaussSeidel(A, b, x0, max_iter, tol)
    tol_reached = 0;
    if nargin < 4
        max_iter = 1e6;
    end
    
    if nargin < 5
        tol = 1e-6;
    end
    x = zeros(max_iter + 1,length(x0));
    x(1,:) = x0;
    lb = length(b);
        for i = 1:max_iter
            for j = 1:lb
                x(i+1,j) = (1/A(j,j))*(b(j) - ( dot(A(j,1:j),x(i+1,1:j)) + dot(A(j,j:lb),x(i,j:lb)) - A(j,j)*( x(i+1,j) + x(i,j) )));
            end

            if norm(x(i+1)-x(i),inf)./norm(x(i),inf) <tol
                tol_reached = 1;
                break;
            end
        end
        x = x(i,:)';
end