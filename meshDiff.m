classdef meshDiff
    %meshDiff Creates the mesh discretization of a flat plate geometry.
    %   Uniform mesh is created in both directions (x-direction and
    %   y-direction) to discretize a flat-plate geometry.

    properties
        material
        n_x
        n_y
        edges_x
        edges_y
        dx
        dy
        x
        y
        A
        Q
        k
    end

    methods
        function obj = meshDiff(material, n_x, n_y)
            %meshDiff Construct an instance of meshDiff class
            %   Basic mesh parameters are defined.
            obj.material = material;
            obj.n_x = n_x; %number of cells in x dimension
            obj.n_y = n_y; %number of cells in y dimension
            obj.edges_x = obj.n_x+1;
            obj.edges_y = obj.n_y+1;
            obj.material = obj.material.calculateProperties(obj.edges_x,obj.edges_y);

        end
        
        function obj = mesh_uniform_sizing(obj)
            %calculating mesh size in both directions
            obj.dx = obj.material.length/(obj.n_x);
            obj.dy = obj.material.height/(obj.n_y);

            %calculating edge cells coordinates - flux is edge-centered,
            %properties are cell-centered
            obj.x = 0:obj.dx:obj.material.length;
            obj.y = 0:obj.dy:obj.material.height;

            %initialize coef matrix A and source vector Q
            obj.A = zeros(obj.edges_x*obj.edges_y, obj.edges_x*obj.edges_y);
            obj.Q = zeros(obj.edges_x*obj.edges_y,1);
            obj.k = zeros(obj.edges_y, obj.edges_x);

            %setting coefficients matrix A for cell that are not boundaries
            for i = 2: obj.edges_y-1
                for j = 2:obj.edges_x-1
                    %giving a unique index for each cell, begins in left
                    %bottom and goes to the right and the to top
                    k = (i-1)*obj.edges_x + j;
                    obj.k(i,j) = k;


                    %applying FVM formulation described in lesson 14
                    %ai-1,j
                    a_left = -(obj.material.D(i,j)*obj.dy + obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                    %ai+1,j
                    a_right = -(obj.material.D(i,j+1)*obj.dy + obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j)*obj.dx + obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j)*obj.dx + obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i,j+1)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                    obj.A(k,k) = a_center;
                    obj.A(k,k-1) = a_left;
                    obj.A(k,k+1) = a_right;
                    obj.A(k,k-obj.edges_x) = a_bottom;
                    obj.A(k,k+obj.edges_x) = a_top;
                    obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4 + obj.material.source(i+1,j)*obj.dx*obj.dy/4 + obj.material.source(i,j+1)*obj.dx*obj.dy/4 + obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;
                end
            end
        end
        function obj = setReflectiveBoundary(obj, edge)
            %setReflectiveBoundary method used to set a boundary as
            %reflective in the edge passed as input to the method.
            if edge == "left"
                index = 1:obj.edges_x:obj.edges_x*obj.edges_y;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k); 

                    %ai-1,j
                    a_left = 0;
                    %ai+1,j
                    a_right = -(obj.material.D(i,j+1)*obj.dy + obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j+1)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    
                    obj.A(k,k) = a_center;
                    obj.A(k,k+1) = a_right;
                    obj.A(k,k-obj.edges_x) = a_bottom;
                    obj.A(k,k+obj.edges_x) = a_top;
                    %add source
                    obj.Q(k) = 0;obj.material.source(i,j+1)*obj.dx*obj.dy/4 + obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;
                end
                %bottom node
                k = index(1);
                [i,j] = obj.k2ij(k); 

                %ai-1,j
                a_left = 0;
                %ai+1,j
                a_right = -(obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    
                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k+obj.edges_x) = a_top;
                obj.Q(k) = obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;

                %top node
                k = index(end);
                [i,j] = obj.k2ij(k); 

                %ai-1,j
                a_left = 0;
                %ai+1,j
                a_right = -(obj.material.D(i,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = -(obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    
                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k-obj.edges_x) = a_bottom;
                obj.Q(k) = obj.material.source(i,j+1)*obj.dx*obj.dy/4;

            elseif edge == "bottom"
                index = 1:1:obj.edges_x;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k);
                    %ai-1,j
                    a_left = -(obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                    %ai+1,j
                    a_right = -(obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = 0;
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j)*obj.dx + obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                    obj.A(k,k) = a_center;
                    obj.A(k,k+1) = a_right;
                    obj.A(k,k-1) = a_left;
                    obj.A(k,k+obj.edges_x) = a_top;
                    %add source
                    obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4 + obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;
                end
                %left node
                k = index(1);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = 0;
                %ai+1,j
                a_right = -(obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;

                %right node
                k = index(end);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -(obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = 0;
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4;
            elseif edge == "top"
                index = obj.edges_x*(obj.edges_y - 1) + 1:1:obj.edges_x*obj.edges_y;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k);
                    %ai-1,j
                    a_left = -(obj.material.D(i,j)*obj.dy)/(2*obj.dx);
                    %ai+1,j
                    a_right = -(obj.material.D(i,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j)*obj.dx + obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = 0;
                    %ai,j
                    a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                    obj.A(k,k) = a_center;
                    obj.A(k,k+1) = a_right;
                    obj.A(k,k-1) = a_left;
                    obj.A(k,k-obj.edges_x) = a_bottom;
                    %add source
                    obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4 + obj.material.source(i,j+1)*obj.dx*obj.dy/4;
                end
                %left node
                k = index(1);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = 0;
                %ai+1,j
                a_right = -(obj.material.D(i,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = -(obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    
                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j+1)*obj.dx*obj.dy/4;

                %right node
                k = index(end);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -(obj.material.D(i,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = 0;
                %ai,j-1
                a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4;
            else
                index =  obj.edges_x:obj.edges_x:obj.edges_x*obj.edges_y;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k); 

                    %ai-1,j
                    a_left = -(obj.material.D(i,j)*obj.dy + obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                    %ai+1,j
                    a_right = 0;
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                    obj.A(k,k) = a_center;
                    obj.A(k,k-1) = a_left;
                    obj.A(k,k-obj.edges_x) = a_bottom;
                    obj.A(k,k+obj.edges_x) = a_top;
                    obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4 + obj.material.source(i+1,j)*obj.dx*obj.dy/4;
                end
                %bottom node
                k = index(1);
                [i,j] = obj.k2ij(k); 

                %ai-1,j
                a_left = -(obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = 0;
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4;

                %top node
                k = index(end);
                [i,j] = obj.k2ij(k); 
                %ai-1,j
                a_left = -(obj.material.D(i,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = 0;
                %ai,j-1
                a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    
                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4;
            end                
        end

        function obj = setVacuumBoundary(obj,edge)
            %setVacuumBoundary method used to set a boundary as
            %vacuum in the edge passed as input to the method.
            if edge == "left"
                index = 1:obj.edges_x:obj.edges_x*obj.edges_y;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k);
                    %ai-1,j
                    a_left = -1/(2.1312*obj.material.D(i,j));
                    %ai+1,j
                    a_right = -(obj.material.D(i,j+1)*obj.dy + obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                    obj.A(k,k) = a_center;
                    obj.A(k,k+1) = a_right;
                    obj.A(k,k-obj.edges_x) = a_bottom;
                    obj.A(k,k+obj.edges_x) = a_top;
                    %add source
                    obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4 + obj.material.source(i+1,j)*obj.dx*obj.dy/4;
                end
                %bottom node
                k = index(1);

                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -1/(2.1312*obj.material.D(i,j));
                %ai+1,j
                a_right = -(obj.material.D(i,j+1)*obj.dy + obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4;

                %top node
                k = index(end);

                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -1/(2.1312*obj.material.D(i,j));
                %ai+1,j
                a_right = -(obj.material.D(i,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4;

            elseif edge == "bottom"
                index = 1:1:obj.edges_x;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k);
                    %ai-1,j
                    a_left = -(obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                    %ai+1,j
                    a_right = -(obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = -1/(2.1312*obj.material.D(i,j));
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j)*obj.dx + obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                    obj.A(k,k) = a_center;
                    obj.A(k,k+1) = a_right;
                    obj.A(k,k-1) = a_left;
                    obj.A(k,k+obj.edges_x) = a_top;
                    %add source
                    obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4 + obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;
                end
                %left node
                k = index(1);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = 0;
                %ai+1,j
                a_right = -(obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = -1/(2.1312*obj.material.D(i,j));
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx + obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;

                %right node
                k = index(end);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -(obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = 0;
                %ai,j-1
                a_bottom = -1/(2.1312*obj.material.D(i,j));
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4;
            elseif edge == "top"
                index = obj.edges_x*(obj.edges_y - 1) + 1:1:obj.edges_x*obj.edges_y;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k);
                    %ai-1,j
                    a_left = -(obj.material.D(i,j)*obj.dy)/(2*obj.dx);
                    %ai+1,j
                    a_right = -(obj.material.D(i,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j)*obj.dx + obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -1/(2.1312*obj.material.D(i,j));
                    %ai,j
                    a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                    obj.A(k,k) = a_center;
                    obj.A(k,k+1) = a_right;
                    obj.A(k,k-1) = a_left;
                    obj.A(k,k-obj.edges_x) = a_bottom;
                    %add source
                    obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4 + obj.material.source(i,j+1)*obj.dx*obj.dy/4;
                
                end
                %left node
                k = index(1);

                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = 0;
                %ai+1,j
                a_right = -(obj.material.D(i,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = -(obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = -1/(2.1312*obj.material.D(i,j));
                %ai,j
                a_center = obj.material.sigma_a(i,j+1)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j+1)*obj.dx*obj.dy/4;
                

                %right node
                k = index(end);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -(obj.material.D(i,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = 0;
                %ai,j-1
                a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = -1/(2.1312*obj.material.D(i,j));
                %ai,j
                a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);
                    

                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4;
            else
                index =  obj.edges_x:obj.edges_x:obj.edges_x*obj.edges_y;
                for a = 2:length(index)-1
                    k = index(a);
                    
                    [i,j] = obj.k2ij(k);
                    %ai-1,j
                    a_left = -(obj.material.D(i,j)*obj.dy + obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                    %ai+1,j
                    a_right = -1/(2.1312*obj.material.D(i,j));
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 + obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                    obj.A(k,k) = a_center;
                    obj.A(k,k-1) = a_left;
                    obj.A(k,k-obj.edges_x) = a_bottom;
                    obj.A(k,k+obj.edges_x) = a_top;
                    %add source
                    obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4 + obj.material.source(i+1,j)*obj.dx*obj.dy/4;
                end
                %bottom node
                k = index(1);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -(obj.material.D(i+1,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = -1/(2.1312*obj.material.D(i,j));
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i+1,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4;

                %top node
                k = index(end);
                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -(obj.material.D(i,j)*obj.dy)/(2*obj.dx);
                %ai+1,j
                a_right = -1/(2.1312*obj.material.D(i,j));
                %ai,j-1
                a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j)*obj.dx*obj.dy/4 - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4;
            end
        end


        function  plotMesh(obj)
            %plotMesh method used to print in command window how the mesh
            %nodes are disposed according to i and j indexes.
            for j = obj.edges_y :-1: 1
                for i = 1 : obj.edges_x
                    fprintf("(%d,%d) = %d ", i,j, (j-1)*(obj.edges_x) + i);
                end
                fprintf("\n")
            end
        end

        function [i,j] = k2ij(obj,k)
            %k2ij function to convert the index node (k) to (i,j) matrix
            %index.
            j = mod(k,obj.edges_x);
            i = floor(k/obj.edges_x)+1;
            if j == 0
                j=obj.edges_x;
                i = i - 1;
            end
        end

        function k = testMesh(obj)
            %testMesh function to test mesh nodes disposal
            k = zeros(obj.edges_y, obj.edges_x);
            for j = 1:obj.edges_x
                for i = 1:obj.edges_y
                    k(i,j) = (i-1)*obj.edges_x + j;
                end
            end
        end
    end
end