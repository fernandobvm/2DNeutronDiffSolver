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
            if nargin == 5 && method == "SOR"
                    obj.phi = SOR(obj.mesh.A, obj.mesh.Q, zeros(length(obj.mesh.Q),1), omega, max_iter, tol);
            else
                max_iter = 1e5;
                tol = 1e-6;
                    if method == "GS"
                        obj.phi = gaussSeidel(obj.mesh.A, obj.mesh.Q, zeros(length(obj.mesh.Q),1), max_iter, tol);
                    else
                        obj.phi = jacobi(obj.mesh.A, obj.mesh.Q, zeros(length(obj.mesh.Q),1), max_iter, tol);
                    end
            end
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end