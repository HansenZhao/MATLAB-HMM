classdef ExtNum
    properties
        p1,
        p2
    end
    
    properties(Dependent)
        value
    end
    
    methods
        function obj = ExtNum(n)
            if ~exist('n','var')
                n = 1;
            end
            if n<0
                obj = ExtNum(-n);
                obj.p1 = -obj.p1;
            else
                tmp = log10(n);
                if isinf(tmp)
                    error('ExtNum init: invalid number');
                end
                tmp2 = floor(tmp);
                tmp1 = tmp - tmp2;
                obj.p1 = 10^tmp1;
                obj.p2 = tmp2;
            end
        end
        
        function v = get.value(obj)
            v = obj.p1*10^obj.p2;
        end
        
        function v = plus(o1,o2)
            factorDiff = o1.p2 - o2.p2;
            v = ExtNum();
            if factorDiff >= 0
                v.p1 = o1.p1 + o2.p1/power(10,factorDiff);
                v.p2 = o1.p2;
            else
                v.p1 = o1.p1/power(10,-factorDiff) + o2.p1;
                v.p2 = o2.p2;
            end
        end
        
        function v = multiValue(obj,n)
            v = ExtNum();
            v.p1 = obj.p1*n;
            v.p2 = obj.p2;
            v = v.reConfigValue();
        end
        
        function v = mtimes(o1,o2)
            v = ExtNum();
            v.p1 = o1.p1*o2.p1;
            v.p2 = o1.p2+o2.p2;
            v.reConfigValue();
        end
        
        function v = mrdivide(o1,o2)
            o2.p1 = 1/o2.p1;
            o2.p2 = -o2.p2;
            v = o1 * o2;
        end
     
        function v = reConfigValue(obj)
            v = ExtNum(obj.p1);
            v.p2 = v.p2 + obj.p2;
        end
        
        function disp(obj)
            fprintf(1,'%.4f x E%d\n',obj.p1,obj.p2);
        end
        
        function b = ge(o1,o2)
            if o1.p2 == o2.p2
                b = o1.p1 >= o2.p1;
            else
                b = o1.p2 >= o2.p2;
            end
        end
    end
    
    methods(Static)
        function v = ExtSum(vec)
            L = length(vec);
            v = vec{1};
            for m = 2:L
                v = v + vec{m};
            end
            v = v.reConfigValue();
        end
        
        function v = cellVecTimes(c1,c2)
            if length(c1)~=length(c2)
                error('ExtNum: cells for vector times should be in the same length, Given %d and %d.',...
                    length(c1),length(c2));
            end
            L = length(c1);
            v = cell(L,1);
            for m = 1:L
                v{m} = c1{m}*c2{m};
            end
            v = ExtNum.ExtSum(v);
        end
        
        function [v,I] = ExtMax(vec)
            v = vec{1};
            I = 1;
            L = length(vec);
            for m = 2:L
                if vec{m} >= v
                    v = vec{m};
                    I = m;
                end
            end
        end
    end
end

