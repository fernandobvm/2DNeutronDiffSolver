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
    end

    methods
        function obj = material(length, height, source, sigma_t)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.length = length;
            obj.height = height;
            obj.source = source;
            obj.sigma_t = sigma_t;
        end


        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end