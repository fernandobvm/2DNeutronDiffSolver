classdef meshDiff
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        materials
        n_x
        n_y
        dx
        dy
        length
        height
    end

    methods
        function obj = meshDiff(material, n_x, n_y)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.materials = material;
            obj.n_x = n_x;
            obj.n_y = n_y;
        end

        function obj = addMaterial(obj,material)
            obj.materials = [obj.materials material];
        end
        
        function obj = generateMeshSBS(obj)
            %Generate mesh for two materials side-by-side with same height
            
            %length
            a = obj.materials(1).length;
            b = obj.materials(2).length;
            va=linspace(0,a,obj.n_x+1);
            vb=linspace(a,b+a,obj.n_x+1);
            obj.length = [va vb(2:end)];

            %height
            c = obj.materials(1).height;
            obj.height = linspace(0, c, obj.n_y+1);

            obj.dx = diff(obj.length);
            obj.dy = diff(obj.height);
        end
    end
end