function [Q,R] = qrfact(A) % from textbook using housholder reflections
    [m, n] = size(A);
    Qt = eye(m);
    for k = 1:n
      z = A(k:m, k);
      w = [ -sign(z(1)) * norm(z) - z(1); -z(2:end) ];
      nrmw = norm(w);
      if nrmw < eps, continue, end
      v = w / nrmw;
      for j = 1:n
        A(k:m, j) = A(k:m, j) - v * (2 * (v' * A(k:m, j)));
      end
      for j = 1:m
        Qt(k:m, j) = Qt(k:m, j) - v * (2 * (v' * Qt(k:m, j)));
      end
    end
    Q = Qt';
    R = triu(A);
end