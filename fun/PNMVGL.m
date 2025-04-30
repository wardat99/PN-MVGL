function [Ak,Ck, Gk,time_cost] = PNMVGL(Xk,P,gamma_1,gamma_2,gamma_3,gamma_4,alpha)
% Perturbation Node-Based Multiview Graph Learning.
% Definition:
%     [Ak,Ck, Gk,time_cost] = CHMVGL(Xk,P,gamma_1,gamma_2,gamma_3,gamma_4,alpha)
%
% Inputs:
% Xk              [n*n*K] cell, contaions the samples for each views.
% gamma_1         scalar, regularization parameter that controls the sparsity of C^k.
% gamma_2         scalar, regularization parameter that controls the degree of C^k.
% gamma_3         scalar, regularization parameter that controls the L_{2,1} norm.
% gamma_4         scalar, regularization parameter that controls the sparsity of C.
% 
% 
% 
% Outputs:
% Ak              3rd-mode tensor [n*n*K], the perturbed learned adjacency matrices.
% Ck              3rd-mode tensor [n*n*K], the perturbed learned Laplacian matrices.
% G               [n*n] matrix, contains the learned perturbed nodes.
% Timecost        scalar, the time cost for each attempt.
%
%

%%%%
%
% Copyright (C)  
%  
%%%%


n=size(Xk{1},1);
[~,K]=size(Xk);
for k=1:K
    Bk(:,:,k)=P'*Xk{k}*Xk{k}'*P;
end

%% Initialization
Ck=zeros(n,n,K);
Vk=Ck; Wk=Ck; Gk=Vk;  Zk=Ck; 
Yk=zeros(n,n,K); Jk=Yk; Mk=Yk; Nk=Ck; Qk=Nk; E=zeros(n-1,n-1); T=zeros(n,n);
I=eye(n); C=T;
max_iter = 2000;
tol=1e-1;


tic
for iter=1:max_iter
    
    %% Update Ek, Zk
    for k=1:K
        Ek(:,:,k)=(2*gamma_1*P'*Zk(:,:,k)*P+alpha*P'*Ck(:,:,k)*P+P'*Yk(:,:,k)*P-Bk(:,:,k)')./(2*gamma_1+alpha);
    end

    for k=1:K
        Uk(:,:,k)=2*gamma_1*P*Ek(:,:,k)*P'+alpha*I.*Ck(:,:,k)-Jk(:,:,k);
        Zk(:,:,k)=(Uk(:,:,k)+(Uk(:,:,k).^2+4*(2*gamma_1+alpha)*gamma_2*I).^0.5)./(4*gamma_1+2*alpha);
        Zk(:,:,k)=I.*Zk(:,:,k);
    end

    %% Update E (\Xi), Ck, Ci (\sigma)
    for k=1:K
        E=(alpha*P'*C*P+P'*T*P)./(2*gamma_4+alpha);
    end

   

    for k=1:K
        Ck_d=I.*((alpha*P*Ek(:,:,k)*P'-Yk(:,:,k)+alpha*Zk(:,:,k)+Jk(:,:,k)+alpha*C+alpha*Vk(:,:,k)+alpha*Wk(:,:,k)-Mk(:,:,k))./(3*alpha));
        Ck_d(Ck_d<0)=0;
        Ck_off=(alpha*P*Ek(:,:,k)*P'-Yk(:,:,k)+alpha*C+alpha*Vk(:,:,k)+alpha*Wk(:,:,k)-Mk(:,:,k))./(2*alpha);
        Ck_off(Ck_off>0)=0;
        Ck(:,:,k)=Ck_d+Ck_off; Ck(:,:,k)=0.5*(Ck(:,:,k)+Ck(:,:,k)');
    end

        
    sum01=0;
    for k=1:K
        sum01=sum01+alpha*Ck(:,:,k)-alpha*Vk(:,:,k)-alpha*Wk(:,:,k)+Mk(:,:,k);
    end


   
          Cii = (sum01+alpha*P*E*P'-T)./((K+1)*alpha);
          Ci_d=Cii.*I; Ci_d(Ci_d<0)=0; 
          Ci_off=Cii-Cii.*I; Ci_off(Ci_off>0)=0; 
          C=Ci_d+Ci_off; C=0.5*(C+C'); clear Cii Ci_d Ci_off sum01
    

    %% Update Wk, Vk, Gk

    for k=1:K
        Wk(:,:,k)=(alpha*Ck(:,:,k)-alpha*C-alpha*Vk(:,:,k)+Mk(:,:,k)+alpha*Vk(:,:,k)'+Nk(:,:,k)')./(2*alpha);
        Vk(:,:,k)=(alpha*Ck(:,:,k)-alpha*C-alpha*Wk(:,:,k)+Mk(:,:,k)+alpha*Wk(:,:,k)'-Nk(:,:,k)+alpha*Gk(:,:,k)+Qk(:,:,k))./(3*alpha);
    end

    for k=1:K
        Gkk(:,:,k)=Vk(:,:,k)-Qk(:,:,k)./alpha;
    end

    [Gk] = l21_operator(Gkk,gamma_3,alpha);
    
        



    %% Update Lag.
    for k=1:K
        dYk(:,:,k)=Ck(:,:,k)-P*Ek(:,:,k)*P';
        Yk(:,:,k)=Yk(:,:,k)+alpha*dYk(:,:,k);

        dJk(:,:,k)=Zk(:,:,k)-I.*Ck(:,:,k);
        Jk(:,:,k)=Jk(:,:,k)+alpha*dJk(:,:,k);

        dMk(:,:,k)=Ck(:,:,k)-C-Vk(:,:,k)-Wk(:,:,k);
        Mk(:,:,k)=Mk(:,:,k)+alpha*dMk(:,:,k);

        dNk(:,:,k)=Vk(:,:,k)-Wk(:,:,k)'; 
        Nk(:,:,k)=Nk(:,:,k)+alpha*dNk(:,:,k);

        dQk(:,:,k)=Gk(:,:,k)-Vk(:,:,k);  
        Qk(:,:,k)=Qk(:,:,k)+alpha*dQk(:,:,k);
    end
    
    
    dT=C-P*E*P'; T=T+alpha*dT;
        
    alpha=1.1*alpha;



    %% stopping criteria
    
    err(1)=norm(dYk(:),'fro'); err(2)=norm(dJk(:),'fro');
    err(3)=norm(dMk(:),'fro'); err(4)=norm(dNk(:),'fro');
    err(5)=norm(dQk(:),'fro'); err(6)=norm(dT,'fro');

    err_max=max(err);
    
    if err_max < tol
       break;
    end
end

time_cost=toc;
%% Extract the perturbed adjacancy matrices from the learned co-hub Laplacian matrices.
ThValueVV=0.003;
        for k=1:K
        Akk=diag(diag(Ck(:,:,k)))-Ck(:,:,k); Akk=Akk./norm(Akk); % perturbed adjacancy matrices
        Akk(Akk>ThValueVV)=1; 
        Akk(Akk<ThValueVV)=0;
        Ak(:,:,k)=Akk; clear Akk
                        
        
        end
        
        ThValue=0.03;
        Aii=diag(diag(C))-C; Aii=Aii/norm(((Aii))); % consensus adjacancy matrix
        Aii(Aii>ThValue)=1; 
        Aii(Aii<ThValue)=0; 
        Ai=Aii; clear Aii

end
