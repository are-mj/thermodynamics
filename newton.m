function x = newton(fun,x0,bounds)
% Solving fun(x) = 0 by newton's method
%  x = newton(fun,x0,[bounds])
%  fun must have the form [f,J] = fun(x), f: function value, J: Jacobian
%    x0 and f must both be n by 1 matrices. J is n by n
%  bounds: optional  n by 2 matrix  
%    bounds(:,1): lower bounds, bounds(:,2) upper bounds
%    Use -Inf and Inf for no bounds for a given element
  
  TolFun = 1e-5;
  TolX   = 1e-4;
  MaxIter = 25;
  x = x0;
  bnd = nargin>2;
  for i = 1:MaxIter
    [f,J] = fun(x);
    dx = J\f(:);
    if bnd
      xnew = x - dx;
      lo = xnew <= bounds(:,1);
      hi = xnew >= bounds(:,2);
      dx(lo) = (x(lo)-bounds(lo,1))/2;
      dx(hi) = (bounds(hi,1) - x(hi))/2;
    end
    x = x - dx;
    if norm(f(:))<TolFun && norm(dx)<TolX
        return
    end
  end
  error(['No convergence after ',num2str(MaxIter),' iteraions'])
end