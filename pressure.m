function p = pressure(th,z)
% Vector of pressures from correspondinng T and v column vectors: z = [T,v]
  T = z(:,1);
  v = z(:,2);
  p = zeros(size(T));
  for i = 1:length(T)
    Tvcalc(th,T(i),v(i));
    p(i) = th.p;
  end
end