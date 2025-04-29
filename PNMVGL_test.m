clear all
clc
close all
warning off;

addpath('fun')

for j=1:50
n=2^7;
p=0.08;
num_views = 2;
pn_pern=0.02;
perturbed_nodes=int8(pn_pern*n);
str = 80/100;
pn_str = int8(n*str);
Adj  = create_ER_Graph(n,p);
[A,rn] = get_perturbed_graph_diff(Adj,perturbed_nodes,num_views,pn_str);
L=diag(sum(Adj))-Adj;
n_signals=700;
noise_amount=0.1;
for v=1:num_views
X{v} = gen_samples_new(A(:,:,v),n_signals,noise_amount,'heat')';
end
[P] = generate_P(n)';
alpha=1;
delta1=4; delta2=30; delta3=10; delta4=4;
gamma_1=delta1/p; gamma_2=delta2*n; gamma_3=delta3/pn_pern; gamma_4=delta4/p;

[Ak,Ck, Gk,time_cost(j)] = PNMVGL(X,P,gamma_1,gamma_2,gamma_3,gamma_4,alpha);

  for k=1:num_views
     [f1v(k,j),~,~] = compute_f(A(:,:,k),Ak(:,:,k));
  end
   f1_view(j) = mean(f1v(:,j));


 end

fprintf('The avg. F-score across views is: %.2f\n' , mean(mean(f1_view)));
fprintf('The avg. TIme cost is: %.2f\n' , mean(mean(time_cost)));
