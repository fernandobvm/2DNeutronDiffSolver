classdef ExamplesLib
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        example
        num
        material
        mesh
        solver
        n_x
        n_y
        X
        Y
        phi
    end

    methods
        function obj = ExamplesLib(num, n_x,n_y)
            %ExamplesLib Run examples provided in Project Instructions File
            %   Each example number represent a combination of both
            %   materials mentioned in the file varying c value accordingly
            obj.num = num;
            obj.n_x = n_x;
            obj.n_y = n_y;

            switch obj.num
                case 1
                    %defining materials
                    mat1 = material(2,2,1,1,0.9);
                    mat2 = material(2,2,0,2,0.5);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);
                case 2
                    %defining materials
                    mat1 = material(2,2,1,1,0.9);
                    mat2 = material(2,2,0,2,0.8);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);
                case 3
                    %defining materials
                    mat1 = material(2,2,1,1,0.9);
                    mat2 = material(2,2,0,2,0.9);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);
                case 4
                    %defining materials
                    mat1 = material(2,2,1,1,0.9);
                    mat2 = material(2,2,0,2,0.99);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);
                case 5
                    %defining materials
                    mat1 = material(2,2,1,1,0.99);
                    mat2 = material(2,2,0,2,0.5);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);
                case 6
                    %defining materials
                    mat1 = material(2,2,1,1,0.99);
                    mat2 = material(2,2,0,2,0.8);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);
                case 7
                    %defining materials
                    mat1 = material(2,2,1,1,0.99);
                    mat2 = material(2,2,0,2,0.9);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);
                case 8
                    %defining materials
                    mat1 = material(2,2,1,1,0.99);
                    mat2 = material(2,2,0,2,0.99);
                    %combining materials
                    obj.material = mat1.addMaterial(mat2);

            end
            %creating the mesh
            obj.mesh = meshDiff(obj.material, obj.n_x, obj.n_y);
            obj.mesh = obj.mesh.mesh_uniform_sizing();
            %setting boundary conditions
            obj.mesh = obj.mesh.setVacuumBoundary("bottom");
            obj.mesh = obj.mesh.setVacuumBoundary("left");
            obj.mesh = obj.mesh.setReflectiveBoundary("right");
            obj.mesh = obj.mesh.setVacuumBoundary("top");
            %creating the solver
            obj.solver = DiffSolver2D(obj.mesh);
        end

        function obj = run(obj)
            %solving matrix equation
            obj.solver = obj.solver.solve();
            %rearranging the data to plot
            [obj.X,obj.Y] = meshgrid(obj.mesh.x, obj.mesh.y);
            obj.phi = reshape(obj.solver.phi,obj.mesh.edges_y,obj.mesh.edges_x)';
            %ploting the results
            figure
            surf(obj.X,obj.Y,obj.phi)
            xlabel("X")
            ylabel("Y")
            zlabel("Neutron Flux")
            title("Example "+num2str(obj.num))
        end

        function obj = compare(obj,num2)
            ex = ExamplesLib(num2,obj.n_x,obj.n_y);
            %solving for first example
            obj.solver = obj.solver.solve();
            %rearranging the data to plot
            [obj.X,obj.Y] = meshgrid(obj.mesh.x, obj.mesh.y);
            obj.phi = reshape(obj.solver.phi,obj.mesh.edges_y,obj.mesh.edges_x)';
            %ploting the results
            figure
            surf(obj.X,obj.Y,obj.phi, 'FaceColor','g', 'FaceAlpha', 0.5)
            hold on
                        
            %solving for second example
            ex.solver = ex.solver.solve();
            ex.phi = reshape(ex.solver.phi,ex.mesh.edges_y,ex.mesh.edges_x)';
            surf(obj.X,obj.Y,ex.phi, 'FaceColor','r', 'FaceAlpha', 0.5)
            xlabel("X")
            ylabel("Y")
            zlabel("Neutron Flux")
            legend("Example "+num2str(obj.num), "Example "+num2str(num2))
            title("Comparing Example "+num2str(obj.num)+ " and Example "+num2str(num2))
        end
    end
end