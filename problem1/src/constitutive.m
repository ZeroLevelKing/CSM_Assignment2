function [Es,Ss,P]=constitutive(B,properties,d,elem,sp)

    [elements,NNE]=size(elem);

    E=properties(1);
    nu=properties(2);
    G=E/2/(1+nu);
    K=2*G/3*(1+nu)/(1-2*nu);
    
    % plane strain
    % lam=2*G/(1-2*nu);
    % D=[lam*(1-nu)     lam*nu         0;
    %     lam*nu        lam*(1-nu)     0;
    %     0                   0        G];

    % plane stress
    lam=E/(1+nu)/(1-nu);
    D=[lam    lam*nu  0;
       lam*nu lam     0;
       0      0       lam*(1-nu)/2];
   
    
    m=[1 1 0];

    Ss=zeros(3*elements,1);
    Es=zeros(3*elements,1);
    P=zeros(elements,1);
    
    for e=1:elements
        b=B{e};
        u=zeros(NNE*sp,1);
        for i=1:NNE
            u(i*sp-1,1)=d(elem(e,i)*sp-1);
            u(i*sp,1)=d(elem(e,i)*sp);
        end
        EE=b*u;
        P(e)=-K*m*EE;

        S=D*EE;

        for i=1:3
            Es(e*3+1-i)=EE(4-i);
            Ss(e*3+1-i)=S(4-i);
        end


    end
    
end