classdef DiffSolver2D
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mesh
        phi
    end

    methods
        function obj = DiffSolver2D(mesh)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.mesh = mesh;
        end

        function obj = solve(obj, method, max_iter, tol, omega)
            if nargin == 1
                max_iter = 1e5;
                tol = 1e-6;
                method = "Jacobi";
            elseif nargin == 2
                max_iter = 1e5;
                tol = 1e-6;
                omega = 1.3;
            elseif nargin == 3
                tol = 1e-6;
                omega = 1.3;
            elseif nargin == 4
                omega = 1.3;
            end
            
            if method == "Jacobi"
                obj.phi = jacobi(sparse(obj.mesh.A), obj.mesh.Q, zeros(length(obj.mesh.Q),1), max_iter, tol);
            elseif method == "GS"
                obj.phi = gaussSeidel(sparse(obj.mesh.A), obj.mesh.Q, zeros(length(obj.mesh.Q),1), max_iter, tol);
            elseif method == "Matlab"
                obj.phi = sparse(obj.mesh.A)\obj.mesh.Q;
            else
                obj.phi = SOR(sparse(obj.mesh.A), obj.mesh.Q, zeros(length(obj.mesh.Q),1), omega, max_iter, tol);
            end
            
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end