function [A,p_n] = get_perturbed_graph_diff(Adj,num_nodes,num_views,str)
% This code is to create a perturbed graph 
% Inputs: Adj >> matrix >> is the adjacacny matrix with dimension n*n
%         num_nodes >> scale >> is the number of views
%
% Outputs: A >> tensor >> is a tensor with dimension (n*n*V) where each 
%          frontal slice is represented as a single perturbed graph 
%          rn >> vector >> is a 1*V vector contains the values of perturbed
%          nodes.




[n,~]=size(Adj);
rng('shuffle');
for v=1:num_views
A(:,:,v)=Adj;
end
rn=[];


for  j=1:num_views 
    x = setdiff(1:n,rn);
    for i=1:num_nodes
        rn(i)=x(randi(numel(x)));
        rv = generate_binary_vector(n, str);
        A(rn(i),:,j)=rv;
        A(:,rn(i),j)=rv';
        p_n{j}=rn;
        x=setdiff(1:n,rn);
    end
end
 



% % for i=1:num_nodes
% %     %rn(i)=x(randi(numel(x)));
% %     for  j=1:num_views 
% %         rn(j)=x(randi(numel(x)));
% %         rv = int32(randi([0, 1], [1, n]));
% %         A(rn(i),:,j)=rv;
% %         A(:,rn(i),j)=rv';
% %         p_n{j}=rn;
% %     end
% %     x=setdiff(1:n,rn);
% % 
% % end



% rn(1)=randi([1,n],1);
% rv1(:,1) = int32(randi([0, 1], [1, n]));
% x = setdiff(1:n, [rn(1)]);
% rn(2) = x(randi(numel(x)));
% rv1(:,2) = int32(randi([0, 1], [1, n]));
% A(:,:,1)=Adj;
% A(rn(1),:,1)=rv1(:,1);
% A(:,rn(1),1)=rv1(:,1)';
% A(rn(2),:,1)=rv1(:,2);
% A(:,rn(2),1)=rv1(:,2)';
% 
% rv2(:,1) = int32(randi([0, 1], [1, n]));
% rv2(:,2) = int32(randi([0, 1], [1, n]));
% A(:,:,2)=Adj;
% A(rn(1),:,2)=rv2(:,1);
% A(:,rn(1),2)=rv2(:,1)';
% A(rn(2),:,2)=rv2(:,2);
% A(:,rn(2),2)=rv2(:,2)';
% 
% % for v=1:num_views
% %         for j=1:size(rn,2)
% %              A(:,:,v)=Adj;
% %              A(rn(j),:,v)=rv(:,j);
% %              A(:,rn(j),v)=rv(:,j)';
% %         end
% % end
% 
% 
% % for v=1:num_views
% %     if v==2
% %         x = setdiff(1:n, [rn(v-1)]);
% %         rn(v) = x(randi(numel(x)));
% %         rv(:,v) = int32(randi([0, 1], [1, n]));        
% %     else
% %         rng('shuffle');
% %         rn(v)=randi([1,n],1);
% %         rv(:,v) = int32(randi([0, 1], [1, n]));  
% %     end
% % end
% % 
% %         for v=1:num_views
% %             A(:,:,v)=Adj;
% %             A(rn(v),:,v)=rv(:,v);
% %             A(:,rn(v),v)=rv(:,v)';
% %         end
function binary_vector = generate_binary_vector(n, str)
    if str > n || str < 0
        error('r must be between 0 and n');
    end
    binary_vector = zeros(1, n);  % Initialize vector with zeros
    indices = randperm(n, str);     % Select str unique positions randomly
    binary_vector(indices) = 1;   % Set selected positions to 1
end

end