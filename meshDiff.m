classdef meshDiff
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        material
        n_x
        n_y
        dx
        dy
        length
        height
        A
        Q
    end

    methods
        function obj = meshDiff(material, n_x, n_y)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.material = material;
            obj.n_x = n_x; %nodes number
            obj.n_y = n_y; %nodes number
        end
        
        function obj = generateMeshSBS(obj)
            %Generate mesh for two materials side-by-side with same height
            
            %length
            a = obj.material(1).length;
            b = obj.material(2).length;
            va=linspace(0,a,obj.n_x+1);
            vb=linspace(a,b+a,obj.n_x+1);
            obj.length = [va vb(2:end)];

            %height
            c = obj.material(1).height;
            obj.height = linspace(0, c, obj.n_y+1);

            obj.dx = diff(obj.length);
            obj.dy = diff(obj.height);
        end

        function obj = mesh_uniform_sizing(obj)
            obj.dx = obj.material.length/(obj.n_x);
            obj.dy = obj.material.height/(obj.n_y);
            obj.length = linspace(0, obj.material.length, obj.n_x);
            obj.height = linspace(0, obj.material.height, obj.n_y);

            sigma_a = obj.material.sigma_t(1)*ones(1,obj.n_x);
            sigma_a(obj.length> obj.length/obj.material.num_materials) = obj.material.sigma_t(2);
            sigma_a(obj.length == obj.length/obj.material.num_materials) = (obj.material.sigma_t(1) + obj.material.sigma_t(2))/2;
            obj.A = zeros((obj.n_x)*(obj.n_y), (obj.n_x)*(obj.n_y));
            obj.Q = zeros(1, (obj.n_x)*(obj.n_y));

            if length(obj.material.source) < 1
                obj.material.source = obj.material.source*ones(obj.n_x+1, obj.n_y+1);
            end
            D = 1;
            for i = 2: obj.n_x-1
                for j = 2:obj.n_y-1
                    k = (j-1)*obj.n_x + i;

                    obj.A(k,k) = 2*D/obj.dx^2 + 2*D/obj.dy^2 + sigma_a(i);
                    obj.A(k,k-1) = -D/obj.dx^2;
                    obj.A(k,k+1) = -D/obj.dx^2;
                    obj.A(k,k-obj.n_x) = -D/obj.dy^2;
                    obj.A(k,k+obj.n_x) = -D/obj.dy^2;
                    obj.Q(k) = 0; %obj.material.source(i,j);
                end
            end
        end



        function obj = plotMesh(obj)
            for j = obj.n_x + 1 :-1: 1
                for i = 1 : obj.n_y+1
                    fprintf("(%d,%d) = %d ", j,i, (j-1)*(obj.n_x+1) + i);
                end
                fprintf("\n")
            end
        end
    end
end