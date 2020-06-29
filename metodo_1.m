warning('off')
pkg load symbolic
function [x_k,itera,ep] = metodo_1(F,vars,x_0,tol) 
    w=[];
    x=[];
    y=[];
    z=[];
    error=[];
    a=1;
    b=-2;
    itera=0;
    iteraciones=[];
    ep=tol+1;
    JF= jacobian(F,vars);

    while(ep>tol)
        Jx_kF = linsolve(double(Fx(JF,vars,x_0)), double(Fx(F, vars, x_0))');

        y_k=x_0 - (1/2)*(Jx_kF');

        JF_sol = linsolve(double(Fx(JF,vars,y_k)), double(-Fx(F, vars, x_0))');

        z_k = x_0 + JF_sol';

        M = double(a * Fx(JF,vars,x_0) + b *  Fx(JF,vars,y_k));

        Fz= double(Fx(F,vars,z_k));

        M_solve = linsolve(M,-Fz');

        x_k = z_k + M_solve';

        x_0=x_k;

        itera = itera+ 1;
        iteraciones=[iteraciones itera];
        w=[w x_0(1)];
        x=[x x_0(2)];
        y=[y x_0(3)];
        z=[z x_0(4)];
        error=[error ep];
        ep=double(norm(Fx(F,vars,x_0)));
    endwhile
    plot(iteraciones,error,'b--o') 
     
endfunction

function F = Fx(f,vars,x_0) 
    F=subs(f, vars, x_0);     
endfunction