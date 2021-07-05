
function [total, count] = patchiness_helper_f_radius(Idx, nearby_points_idr_idx, X, pcutoff, ncutoff)

count=0; %counts the number of interactions analyzed
total=0; %totals the correlation param.


for i=1:size(Idx,1)
    relevant_idx = Idx{i}; 
    potential_i = X(i,4); 
    if ~isempty(Idx{i})
        Idx{i} = sort(Idx{i});
        for j = 1:length(relevant_idx)
            idx = relevant_idx(j); 
            mapped_idx = nearby_points_idr_idx(idx);
            %disp(idx);
            %disp(i);
            %if mapped_idx ~= i && (((i > mapped_idx) && isempty(Idx{mapped_idx})) || i < mapped_idx)
            if mapped_idx ~= i %we are double counting but thats fine
                count=count+1; 
                potential_j = X(mapped_idx,4); 
                if (potential_i >=pcutoff && potential_j >= pcutoff) ||(potential_i <=ncutoff && potential_j <= ncutoff)
                    total=total+1;
                    %adds 1 to the correlation total if the neighbors are of the same charge
                elseif (potential_i >=pcutoff && potential_j <= ncutoff) ||(potential_i <=ncutoff && potential_j >=pcutoff)
                    total=total-1;
                    %subtracts 1 from the total if the neighbors are of the opposite charge
                else
                    total=total+0;
                end
            end
        end
    end 
end

end