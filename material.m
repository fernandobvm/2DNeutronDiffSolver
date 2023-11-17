classdef material
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        length
        height
        source
        sigma_t
        geometry
        num_materials;
        D
        sigma_s
        sigma_a
    end

    methods
        function obj = material(length, height, source, sigma_t, c)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.length = length;
            obj.height = height;
            obj.source = source;
            obj.sigma_t = sigma_t;
            obj.sigma_s = c*sigma_t;
            obj.D = 1/(3*sigma_t);% isotropic scattering
            obj.sigma_a = obj.sigma_t-obj.sigma_s;
            obj.num_materials = 1;
        end

        function obj = addMaterial(obj, material)
            %addMaterial Add a different material on the side with same
            %length and height of first material
            obj.length = 2*obj.length;
            obj.height = obj.height;
            obj.source = [obj.source material.source];
            obj.D = [obj.D material.D];
            obj.sigma_a = [obj.sigma_a material.sigma_a];
            obj.num_materials = obj.num_materials +1;

        end

        function obj = calculateProperties(obj,edges_x,edges_y)
            if obj.num_materials == 2
                D = ones(edges_y, edges_x);
                sigma_a = ones(edges_y, edges_x);
                source = ones(edges_y, edges_x);

                if mod(edges_x,2) == 0
                    D(:,1:edges_x/2) = obj.D(1);
                    D(:,edges_x/2+1:end) = obj.D(2);
                    sigma_a(:,1:edges_x/2) = obj.sigma_a(1);
                    sigma_a(:,edges_x/2+1:end) = obj.sigma_a(2);
                    source(:,1:edges_x/2) = obj.source(1);
                    source(:,edges_x/2+1:end) = obj.source(2);

                    obj.D = D;
                    obj.sigma_a = sigma_a;
                    obj.source = source;
                else
                    D(:,1:floor(edges_x/2)) = obj.D(1);
                    D(:,ceil(edges_x/2)+1:end) = obj.D(2);
                    D(:,ceil(edges_x/2)) = (obj.D(1)+obj.D(2))/2;
                    sigma_a(:,1:floor(edges_x/2)) = obj.sigma_a(1);
                    sigma_a(:,ceil(edges_x/2)+1:end) = obj.sigma_a(2);
                    sigma_a(:,ceil(edges_x/2)) = (obj.sigma_a(1)+obj.sigma_a(2))/2;
                    source(:,1:floor(edges_x/2)) = obj.source(1);
                    source(:,ceil(edges_x/2)+1:end) = obj.source(2);
                    source(:,ceil(edges_x/2)) = (obj.source(1)+obj.source(2))/2;
                    
                    obj.D = D;
                    obj.sigma_a = sigma_a;
                    obj.source = source;
                end
            else
                obj.D = obj.D*ones(edges_y, edges_x);
                obj.sigma_a = obj.sigma_a*ones(edges_y, edges_x);
                obj.source = obj.source*ones(edges_y, edges_x);
            
            end
            
        end


        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end