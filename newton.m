function x = newton(fun,x0)
% Solving fun(x) = 0 by Newton's method
% fun: [f,J] = fun(x), where f: function value, J: Jacobian
% x0 : Start point for search
%
% More robust alternative: fsolve from the optimization toolbox
  
  TolFun = 1e-5;
  TolX   = 1e-4;
  MaxIter = 25;
  x = x0;
  for i = 1:MaxIter
    [f,J] = fun(x);
    dx = J\f(:);
    x = x - dx;
    if norm(f(:))<TolFun && norm(dx)<TolX
        return
    end
  end
  error(['No convergence after ',num2str(MaxIter),' iteraions'])
end