warning('off')
pkg load symbolic
function [x_k, itera, ep] = metodo_2(F, vars, x_0, tol)
  w = [];
  x = [];
  y = [];
  z = [];
  error = [];
  itera = 0;
  iteraciones = [];
  ep = tol + 1;
  
  JF = jacobian(F, vars);
  
  while(ep > tol)

    mdx = linsolve(double((Fx(JF, vars, x_0))), double(Fx(F, vars, x_0))');
    
    y_k = x_0 - (1/2)*(mdx');
    
    z_k = (1/3)*(4*y_k - x_0);
  
    mxz = double((Fx(JF, vars, x_0)) - 3 * (Fx(JF, vars, z_k)));
 
    mxzxSolve = linsolve(mxz, double(Fx(F, vars, x_0))');
    
    u_k = y_k + mxzxSolve';
  
    mxz2 = double((Fx(JF, vars, x_0) - 3 * Fx(JF, vars, z_k)));
    
    
    mxzuSolve = linsolve(2*mxz2, double(Fx(F, vars, x_0))');
    
    v_k = u_k + mxzuSolve';
  
    mxzvSolve = linsolve(2*mxz2, double(Fx(F, vars, v_k))');
    
    x_k = v_k + mxzvSolve';
  
    x_0 = x_k;
  
    itera = itera + 1;
    iteraciones = [iteraciones itera];
    w = [w x_0(1)];
    x = [x x_0(2)];
    y = [y x_0(3)];
    z = [z x_0(4)];
    error = [error ep];
    ep = double(norm(Fx( F, vars, x_0)));
  endwhile
  plot(iteraciones, error, 'b--o');
 
endfunction

function F = Fx(f, vars, x_0)
  F = subs(f, vars, x_0);
endfunction
