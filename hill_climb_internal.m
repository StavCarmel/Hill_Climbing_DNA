function [dna,obj_score,iter_num] = hill_climb_internal(AAseq)
% this function will extract the first and final sequences and scores in order to
% build the table of the final answer easily.
% in order to extract only the max seq and max score use the 'hill_climb" function
load('aa_codons_map.mat');

% Initiation:
iter_num=0;
changed=ones(1,length(AAseq)-1); % the position of an aa will become 0 if it had changed and improved the sequence
% -1: because can't change the start codon
iteration_dna=aa2nt(AAseq); %gives the first guess randomly
iteration_score=dracarys(iteration_dna);

while any(changed) % if there is still a position that contains 1 it means that it still didn't change to improve the sequnce and the algorithm can proceed 
      iter_num=iter_num+1;
      % selecting a random amino acid
      aa_change=randi(length(AAseq)); %get the position of the aa that we will change randomly
      if aa_change~=1 %can?t change the start codon
          % calculating the score for all the neighbours
          if changed(aa_change-1)==1 %if it is 1 it means that this aa can still change and improve
             codons_options=values(aa_codons_map,{AAseq(aa_change)});
             dna_neighbors=cell(length(codons_options{1})-1,1);
             score_neighbor=zeros(length(codons_options{1})-1,1);
             count=0;
             for codon=1:length(codons_options{1})
                 cur_dna=iteration_dna(end,:);
                 cur_codon=cur_dna((aa_change*3)-2:aa_change*3);
                 if ~strcmp(cur_codon,codons_options{1}{codon}) 
                    count=count+1;
                    dna_neighbors{count,1}=cur_dna;
                    dna_neighbors{count,1}((aa_change*3)-2:aa_change*3)=string(codons_options{1}{codon});
                    score=dracarys(dna_neighbors{count,1});
                    score_neighbor(count,1)=score;
                 end
             end
          end
          [max_neighbor,ind]=max(score_neighbor);
           if max_neighbor<=iteration_score(end) %check if the max from this itteration is higher than the last max
              changed(aa_change-1)=0; % change to 0 because we didn't succed to improve the score with this position
              continue
           end
           iteration_dna(end+1,:)=dna_neighbors{ind,:};
           iteration_score(end+1,1)=max_neighbor;
           changed=ones(1,length(AAseq)-1);
      end
end
dna={iteration_dna(1,:);iteration_dna(end,:)};
obj_score=[iteration_score(1);iteration_score(end)];
end




