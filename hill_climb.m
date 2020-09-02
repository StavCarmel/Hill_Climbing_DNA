function [max_dna,max_obj_score,iter_num] = hill_climb(AAseq)
[dna,obj_score,iter_num] = hill_climb_internal(AAseq);
max_dna=dna(end,:);
max_obj_score=obj_score(end,1);

end

