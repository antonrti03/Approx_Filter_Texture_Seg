function [REC, SM, NS]=Interp_local(I, Ey, Ex, Ez, Ei)
% Created by Tuan Nguyen, 02/04/2021
%
[Y,X] = size(I);
SM = zeros(Y,X);
REC = zeros(Y,X);
SY=zeros(1,Y*X); SX=zeros(1,Y*X); 
y=0;
while (y<Y)
      x=0;      
      while (x<X)
          if SM(y+1,x+1)==0
              NS=I(y+1,x+1); SN=abs(NS);
              options   = rbfcreate([Ex(SN, 1:Ei(SN)-1); Ey(SN, 1:Ei(SN)-1)], Ez(SN, 1:Ei(SN)-1),'RBFFunction', 'invquadratic', 'RBFConstant', 8);
              phi       = options.('rbfphi');
              rbfconst  = options.('RBFConstant');
              nodes     = options.('x');
              rbfcoeff  = (options.('rbfcoeff'))';
              [dim n]   = size(nodes);
              r         = zeros(1, n);
              %
              s = 0;
              r =  ([x+1; y+1]*ones(1,n)) - nodes;%%%[cols; rows]
              r = sqrt(sum(r.*r, 1));
              s = sum(rbfcoeff(1:n).*feval(phi, r, rbfconst))+rbfcoeff(n+1)+rbfcoeff(n+2)*(x+1)+rbfcoeff(n+3)*(y+1);
              REC(y+1,x+1)=s;
              %
              SM(y+1,x+1)=NS; 
              SP=1;
              SY(SP)=y;
              SX(SP)=x; 
             while (SP>0)                         
                   y1=SY(SP);
                   x1=SX(SP);
                   SP=SP-1;
                for j=-1:1
                    for i=-1:1
                        if (y1+j+1>0) && (x1+i+1>0) && (y1+j+1<=Y) && (x1+i+1<=X) && ...
                            SM(y1+j+1,x1+i+1)==0 && I(y1+1,x1+1)==I(y1+j+1,x1+i+1)                     
                            SM(y1+j+1,x1+i+1)=NS;
                            SP=SP+1;
                            SY(SP)=y1+j;
                            SX(SP)=x1+i;
                            %
                            s = 0;
                            r =  ([x1+i+1; y1+j+1]*ones(1,n)) - nodes;%%%[cols; rows]
                            r = sqrt(sum(r.*r, 1));
                            s = sum(rbfcoeff(1:n).*feval(phi, r, rbfconst))+rbfcoeff(n+1)+rbfcoeff(n+2)*(x1+i+1)+rbfcoeff(n+3)*(y1+j+1);
                            REC(y1+j+1,x1+i+1)=s;
                            %
                        end                                     
                    end
                end                               
             end
          end 
     x=x+1;
      end       
y=y+1;
end
%
end