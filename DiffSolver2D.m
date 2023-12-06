classdef DiffSolver2D
    %DiffSolver2D Solver for 2D Neutron Diffusion Equation
    %   The class receives the A and b matrix obtained from Finite Volumes
    %   Method.

    properties
        mesh
        phi
        converged
    end

    methods
        function obj = DiffSolver2D(mesh)
            %DiffSolver2D Solve the matrix system Ax=b using differente
            %methods.
            %   Jacobi, Gauss-Seidel and SOR methods can be used to solve
            %   the system. If "Matlab" is passed as input for solve method
            %   the system is solved using \ matlab function.
            obj.mesh = mesh;
        end

        function obj = solve(obj, method, max_iter, tol, omega)
            %Method use to solve the equations. max_iter is the parameter
            %to set the maximum number of iterations; tol is the parameter
            %to set the tolerance level; omega is the parameter to set the
            %omega value for SOR method.
            if nargin == 1
                max_iter = 1e6;
                tol = 1e-6;
                method = "Jacobi";
            elseif nargin == 2
                max_iter = 1e6;
                tol = 1e-6;
                omega = 1.3;
            elseif nargin == 3
                tol = 1e-6;
                omega = 1.3;
            elseif nargin == 4
                omega = 1.3;
            end
            
            if method == "Jacobi"
                [obj.phi, obj.converged] = jacobi(sparse(obj.mesh.A), sparse(obj.mesh.Q), zeros(length(obj.mesh.Q),1), max_iter, tol);
            elseif method == "GS"
                [obj.phi, obj.converged] = gaussSeidel(sparse(obj.mesh.A), sparse(obj.mesh.Q), zeros(length(obj.mesh.Q),1), max_iter, tol);
            elseif method == "GS2"
                [obj.phi, obj.converged] = gauss_seidel2(sparse(obj.mesh.A), obj.mesh.Q, zeros(length(obj.mesh.Q),1), max_iter, tol);
            elseif method == "Matlab"
                obj.phi = sparse(obj.mesh.A)\obj.mesh.Q;
            else
                [obj.phi, obj.converged] = SOR(sparse(obj.mesh.A), sparse(obj.mesh.Q), zeros(length(obj.mesh.Q),1), omega, max_iter, tol);
            end
            
        end
    end
end