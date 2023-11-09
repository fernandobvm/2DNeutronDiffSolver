function x = jacobi(A, b, x0, max_iter, tol)
    if nargin < 4
        max_iter = 1000;
    end
    
    if nargin < 5
        tol = 1e-6;
    end
    for i = 1:max_iter
        x = (1./diag(A)).*(b-A*x0+diag(A).*x0);
        if norm(x-x0,inf) < tol
            x0 = x;
            break;
        end
        x0 = x;
    end


end