% a simple 2D circular queue structure
% replaces the oldest element when a new element is added to a full queue
classdef circular_queue
    
    properties
        % parameters
        capacity = [];
        size = [];
        front = [];
        array = [];
    end
    
    methods
        
        % constructor must be able to handle the no input case for object array pre-allocation
        function obj = circular_queue(capacity)
            obj.size = 0;
            obj.front = 0;
            if nargin > 0
                obj.capacity = capacity;
                obj.array = zeros(capacity,2);
            end
        end
        
        function [obj] = push_back(obj,ele)
            obj.front = mod(obj.front+1,obj.capacity+1);
            if obj.front == 0
                obj.front = 1;
            end
            
            % fprintf(1,'front: %d\n',obj.front); % DEBUG
                        
            obj.array(obj.front,:) = ele(:);
            
            if obj.size < obj.capacity
                obj.size = obj.size + 1;
            end
            
            % obj.array % DEBUG
        end
        
        function [obj,val] = mean(obj)
            val = zeros(1,2);
            if obj.size ~= 0
                for i = 1:obj.size
                    val = val + obj.array(i,:);
                end
                val = val / obj.size;
            end
            % val % DEBUG
        end
        
    end % methods
end % class

