classdef meshDiff
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

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
    end

    methods
        function obj = meshDiff(material, n_x, n_y)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.material = material;
            obj.n_x = n_x; %number of cells in x dimension
            obj.n_y = n_y; %number of cells in y dimension
            obj.edges_x = obj.n_x+1;
            obj.edges_y = obj.n_y+1;
            obj.material = obj.material.calculateProperties(obj.edges_x,obj.edges_y);

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
            %considering 1 material only
            %calculating mesh size in both directions
            obj.dx = obj.material.length/(obj.n_x);
            obj.dy = obj.material.height/(obj.n_y);
            % %calculating center cells coordinates
            % obj.x = 0+obj.dx/2:obj.dx:obj.material.length-obj.dx/2;
            % obj.y = 0+obj.dy/2:obj.dy:obj.material.length-obj.dy/2;

            %calculating edge cells coordinates - flux is edge-centered,
            %properties are cell-centered
            obj.x = 0:obj.dx:obj.material.length;
            obj.y = 0:obj.dy:obj.material.height;

            % %calculating D and Sigma_a as constant matrix according cells
            % %coordinates
            % obj.material.D = obj.material.D*ones(obj.edges_y, obj.edges_x);
            % obj.material.sigma_a = obj.material.sigma_a*ones(obj.edges_y, obj.edges_x);
            % 
            % %calculating Source matrix as constant
            % obj.material.source = obj.material.source*ones(obj.edges_y, obj.edges_x);
            % 
            
            %initialize coef matrix A and source vector Q
            obj.A = zeros(obj.edges_x*obj.edges_y, obj.edges_x*obj.edges_y);
            size(obj.A)
            obj.Q = zeros(obj.edges_x*obj.edges_y,1);

            %setting coefficients matrix A for cell that are not boundaries
            for j = 2: obj.edges_x-1
                for i = 2:obj.edges_y-1
                    %giving a unique index for each cell, begins in left
                    %bottom and goes to the right and the to top
                    k = (i-1)*obj.edges_x + j;


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
                    obj.A(k,k-obj.n_x) = a_bottom;
                    obj.A(k,k+obj.n_x) = a_top;
                    %add source
                    obj.Q(k) = obj.material.source(i,j+1)*obj.dx*obj.dy/4 + obj.material.source(i+1,j+1)*obj.dx*obj.dy/4;
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
                obj.A(k,k+obj.n_x) = a_top;
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
                obj.A(k,k-obj.n_x) = a_bottom;
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
                    obj.A(k,k+obj.n_x) = a_top;
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
                obj.A(k,k+obj.n_x) = a_top;
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
                obj.A(k,k+obj.n_x) = a_top;
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
                    obj.A(k,k-obj.n_x) = a_bottom;
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
                obj.A(k,k-obj.n_x) = a_bottom;
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
                obj.A(k,k-obj.n_x) = a_bottom;
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
                    obj.A(k,k-obj.n_x) = a_bottom;
                    obj.A(k,k+obj.n_x) = a_top;
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
                obj.A(k,k+obj.n_x) = a_top;
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
                obj.A(k,k-obj.n_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4;
            end                
        end

        function obj = setVacuumBoundary(obj,edge)
            if edge == "left"
                index = 1:obj.edges_x:obj.edges_x*obj.edges_y;
                for a = 2:length(index)-1
                    k = index(a);
                    [i,j] = obj.k2ij(k);
                    %ai-1,j
                    a_left = -1/(2*2.1312*obj.material.D(i,j));
                    %ai+1,j
                    a_right = -(obj.material.D(i,j+1)*obj.dy + obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

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
                a_left = -1/(2*2.1312*obj.material.D(i,j));
                %ai+1,j
                a_right = -(obj.material.D(i,j+1)*obj.dy + obj.material.D(i+1,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k+1) = a_right;
                obj.A(k,k+obj.edges_x) = a_top;
                %add source
                obj.Q(k) = obj.material.source(i+1,j)*obj.dx*obj.dy/4;

                %top node
                k = index(end);

                [i,j] = obj.k2ij(k);
                %ai-1,j
                a_left = -1/(2*2.1312*obj.material.D(i,j));
                %ai+1,j
                a_right = -(obj.material.D(i,j+1)*obj.dy)/(2*obj.dx);
                %ai,j-1
                a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

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
                    a_bottom = -1/(2*2.1312*obj.material.D(i,j));
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j)*obj.dx + obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

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
                a_bottom = -1/(2*2.1312*obj.material.D(i,j));
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx + obj.material.D(i+1,j+1)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

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
                a_bottom = -1/(2*2.1312*obj.material.D(i,j));
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

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
                    a_top = -1/(2*2.1312*obj.material.D(i,j));
                    %ai,j
                    a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);
                    

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
                a_top = -1/(2*2.1312*obj.material.D(i,j));
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);
                    

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
                a_top = -1/(2*2.1312*obj.material.D(i,j));
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);
                    

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
                    a_right = -1/(2*2.1312*obj.material.D(i,j));
                    %ai,j-1
                    a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                    %ai,j+1
                    a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                    %ai,j
                    a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

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
                a_right = -1/(2*2.1312*obj.material.D(i,j));
                %ai,j-1
                a_bottom = 0;
                %ai,j+1
                a_top = -(obj.material.D(i+1,j)*obj.dx)/(2*obj.dy);
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

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
                a_right = -1/(2*2.1312*obj.material.D(i,j));
                %ai,j-1
                a_bottom = -(obj.material.D(i,j)*obj.dx)/(2*obj.dy);
                %ai,j+1
                a_top = 0;
                %ai,j
                a_center = obj.material.sigma_a(i,j) - (a_left + a_right + a_bottom + a_top);

                obj.A(k,k) = a_center;
                obj.A(k,k-1) = a_left;
                obj.A(k,k-obj.edges_x) = a_bottom;
                %add source
                obj.Q(k) = obj.material.source(i,j)*obj.dx*obj.dy/4;
            end
        end


        function  plotMesh(obj)
            for j = obj.edges_y :-1: 1
                for i = 1 : obj.edges_x
                    fprintf("(%d,%d) = %d ", i,j, (j-1)*(obj.edges_x) + i);
                end
                fprintf("\n")
            end
        end

        function [i,j] = k2ij(obj,k)
            j = mod(k,obj.edges_x);
            i = floor(k/obj.edges_x)+1;
            if j == 0
                j=obj.edges_x;
                i = i - 1;
            end
        end
    end
end