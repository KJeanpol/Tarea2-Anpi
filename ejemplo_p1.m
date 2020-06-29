pkg load symbolic
pkg load dataframe
syms x y z w
f=[(w**2)+x-(3*y)+(4*z)+(3/4),(3*w**2)+x-(y**2)+(z**2)+(13/4),(5*w)+(3*x**2)+y-(4*z**2)-(99/2),(8*w**2)-(14*x)+(6*y**2)-(7*z**2)+7];
vars=[w,x,y,z];
x0=[2 1 2 1];
tol=0.1;
[x_k2, itera2,ep2] = metodo_2(f, vars, x0, tol)
[x_k1,itera1,ep1] = metodo_1(f,vars,x0,tol) 
[x_kNR,iteraNR,epNR] = newton_raphson_nl(f,vars,x0,tol,1000)
C = {"Metodo", "ValorI", "X1","X2","X3","X4", "Iteraciones","Error";  
     "Metodo1", "[2 1 2 1]", x_k1(1), x_k1(2), x_k1(3), x_k1(4), itera1,ep1;
     "Metodo2", "[2 1 2 1]", x_k2(1), x_k2(2), x_k2(3), x_k2(4), itera2,ep2;
     "MetodoNR", "[2 1 2 1]", x_kNR(1), x_kNR(2), x_kNR(3), x_kNR(4), iteraNR,epNR;
};
    dataframe (C)