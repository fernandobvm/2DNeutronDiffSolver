function [x,tol_reached] = jacobi(A, b, x0, max_iter, tol)
    tol_reached = 0;
    if nargin < 4
        max_iter = 1e6;
    end
    
    if nargin < 5
        tol = 1e-6;
    end
    for i = 1:max_iter
        x = (1./diag(A)).*(b-A*x0+diag(A).*x0);
        if norm(x-x0,inf)./x < tol
            x0 = x;
            tol_reached = 1;
            break;
        end
        x0 = x;
    end


end