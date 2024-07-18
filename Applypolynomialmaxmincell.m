function result = Applypolynomialmaxmincell(p, A)
    
    d=size(p,2);
    intermediate_results = cell(1, d);
    intermediate_results{1}=min(maxminId(size(A,1)),p(2,1));

    for i = 2:d
        intermediate_results{i}=min(MaxMinpMatPower(A,p(1,i)),p(2,i));
    end

    result = intermediate_results{1};
    for i = 2:length(intermediate_results)
    result = max(result, intermediate_results{i});
    end

end