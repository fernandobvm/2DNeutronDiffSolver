classdef material
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        length
        height
        source
        sigma_f
        sigma_t
        geometry
        num_materials;
    end

    methods
        function obj = material(length, height, source, sigma_t)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.length = length;
            obj.height = height;
            obj.source = source;
            obj.sigma_t = sigma_t;
            obj.num_materials = 1;
        end


        function obj = addMaterial(obj, source, sigma_t)
            %addMaterial Add a different material on the side with same
            %length and height of first material
            obj.length = (obj.num_materials + 1)*obj.length/obj.num_materials;
            obj.num_materials = obj.num_materials + 1;
            obj.source = [obj.source source];
            obj.sigma_t = [obj.sigma_t sigma_t];
        end


        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end