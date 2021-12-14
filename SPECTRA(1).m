
clear all                                         

ev0 = 8.854187817*1e-12;                
uv0 = 4*pi*1e-7;                        
h=6.6260693*1e-34;
c   = 1/sqrt(ev0*uv0);                  
f0  = [300:1:1000]*1e12;              
degree=0.00;                            
m=0.4*9.11*1e-31;                  
n=6*1e26;                              
p=20*1e-4;                             
e0=1.6*1e-19;                          
wpITO=sqrt(n*e0*e0/m/ev0);             
gamaITO=e0/p/m;   
w0     = 2*pi*f0;                     
lambda = c./f0;
theta  = degree*pi/180;                 
F      = f0/1e12;
W      = w0/1e12;
Lambda = lambda*1e9;


layer='AB';

direction='F';                         
if direction=='F'
    layer=layer;
end
if direction=='B'
    layer=layer(length(layer):-1:1);
end

type_wave='TE';


for kw=1:length(w0)
    
    w=w0(kw);
    k0=w/c;
     para_A=[1^2   1^2 1*1e-9];
     para_B=[1^2  1^2 1*1e-9];
     para_C=[1^2  1^2 1*1e-9];
     para_E=[1^2 1^2 1*1e-9];
     para_D=[1^2  1^2 1*1e-9];

    er0=1.0; ur0=1;                        
    ern=1.0; urn=1;                        

    para=[ para_A;para_B;para_C;para_D;para_E];
    para_er    = para(:,1);             
    para_ur    = para(:,2);             
    para_thick = para(:,3);            

    for k=1:length(layer)
        type=uint8(layer(k))-64;
        er(k)    = para_er(type);
        ur(k)    = para_ur(type);
        thick(k) = para_thick(type);
    end

    C_theta0=sqrt(1-sin(theta).^2);
    C_thetan=sqrt(1-sin(theta).^2*(er0*ur0)/(ern*urn));
    if type_wave=='TE'                  
        pr0=sqrt(er0)/sqrt(ur0)*C_theta0;
        prn=sqrt(ern)/sqrt(urn)*C_thetan;
    end
    if type_wave=='TM'                 
        pr0=sqrt(er0)/sqrt(ur0)/C_theta0;
        prn=sqrt(ern)/sqrt(urn)/C_thetan;
    end

    M=eye(2);
    for k=1:length(layer)
        C_theta=sqrt(1-sin(theta).^2*(er0*ur0)/(er(k)*ur(k)));
        b=(w/c)*thick(k)*sqrt(er(k))*sqrt(ur(k))*C_theta;
        if type_wave=='TE'              
            pr=sqrt(er(k))/sqrt(ur(k))*C_theta;
        end
        if type_wave=='TM'             
            pr=sqrt(er(k))/sqrt(ur(k))/C_theta;
        end
        M1=[cos(b) -i/pr.*sin(b); -i*pr.*sin(b) cos(b)];
        M=M*M1;                   
    end

    m11=M(1,1);m12=M(1,2);m21=M(2,1);m22=M(2,2);
    Eout  = sqrt(1);                                   
    Ein   = ((m11+m12*prn)+(m21+m22*prn)/pr0)/2*Eout;
    Erl   = ((m11+m12*prn)-(m21+m22*prn)/pr0)/2*Eout;
    r(kw) = Erl/Ein;                  
    t(kw) = sqrt(prn/pr0)*Eout/Ein;     
  
end

R  = abs(r).^2;                         
T  = abs(t).^2;                         
A  = 1-R-T;                             



plot(Lambda,R,'b')
title('Reflection')
legend('Reflection ')
xlabel('\lambda(nm)'),ylabel('Reflection')%,'Reflection','Absorption');
axis([200 900 0 1])
hold on