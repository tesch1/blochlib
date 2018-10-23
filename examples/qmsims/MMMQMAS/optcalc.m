function [A,M,S,F,B,J,K,I,C,original] = optcalc(m,n);
% function [A,M,S,F,B,J,K,I,C,original] = optcalc(m,n);
%
% calculates most efficient method of creating the propagators
%
% input:
%  m = number of subpropagators (subP) per cycle [t1 proppagator]
%  n = number of subPs in a t2 propagator
%
% output:
%  A = cell array with formulas for the optimal calculations (doesn't
%     include the calculation of a full cycle (C1))
%  M = number of matrix multiplications (MM) [initial,optimal]
%  S = number of extra matrix assignments for optimal calculation
%  F = vector of partial forward sequences needed
%  B = vector of partial backward sequences needed
%  J = vector of partial forward inverse subP sequences needed
%  K = vector of partial backward inverse subP sequences needed
%  I = vector of individual inverse subP needed
%  C = vector of needed complete cycles
%
% for F,B,J,K: 0 = not needed, otherwise = row of A where used
% for C: [1,N-2,N-1,N], 0 for n/a, otherwise N-2, N-1, or N
%
% labels definitions (character, number, and description):
%  U = 0 = the individual subPs
%  I = 1 = inverse subP: Ix = ctranspose(Ux);
%  F = 2 = partial forward sequences: F3 = U3*U2*U1;
%  B = 3 = partial backwards sequences: R5 = Um*Um-1*Um-2*...*U5;
%  J = 4 = partial forwards inverse sequences: J3 = I1*I2*I3;
%  K = 5 = partial backwards inverse sequences: K5 = I5*I6*I7*...*Im;
%  C = 6 = complete cycles: C1 = Um*Um-1*Um-2*...*U1; C2 = C1*C1; C3 = C1*C2;
%  Z = 7 = element of A cell array
%

% author:  John W. Logan
% version: 19 February 2003

% check inputs...
if gcd(m,n) ~= 1
   error('gcd(m,n) must equal 1')
end
if prod(size(m)) ~= 1
   error('m must be an integer')
end
if prod(size(n)) ~= 1
   error('n must be an integer')
end
if m < 2
   error('m must be greater than 1')
end
if n < 2
   error('n must be greater than 1')
end
if m > floor(2^51)
   error('m must be less than 450359962737049')
end
if n > floor(2^51)
   error('n must be less than 450359962737049')
end

offset = 10^(ceil(log10(m+1))); % label offset

% initialize variables output variables (except A)
M = zeros([1,2]);
S = 0;
F = zeros([1,m]);
B = zeros([1,m]);
J = zeros([1,m]);
K = zeros([1,m]);
I = zeros([1,m]);
C = [1, 0, 0]; % [C1, C(N-1), CN]

% original / initialization
h = 1;
A = cell([m,1]); % t2prop(j)
for j = 1:m
   A{j} = zeros([1,n]);
   for k = n:-1:1
      A{j}(k) = (mod(h-1,m) + 1);
      h = h + 1;
   end
end
original = A;
M(1) = m*n - 1; % (m-1 [full cycle]; m*(n-1) [rational fractions])

% Optimization
%  4 cases:
%   m = 2;
%   n+1 = k*m, k = 2,3,4,...
%   n-1 = k*m, k = 2,3,4,...
%   everything else
if m == 2
   % (offset = 10)
   N = (n-1)/m;
   A{1} = [1, repmat(61,[1,N])];
   A{2} = [11, 71, 2];
   M(2) = N+3; % (AA:1; A:N+2)
elseif (~mod(n+1,m) & (n+1 >= 2*m))
   N = (n+1)/m;
   A{1} = [m, ones(1,N)] + [1, repmat(6,[1,N])]*offset;
   for j = 2:m
      A{j} = [m+1-j, j-1, m+2-j] + [1, 7, 0]*offset;
   end
   M(2) = 3*(m-1)+N; % (AA:m-1; A:N+2*(m-1))
   I = ones([1,m]);
elseif (~mod(n-1,m) & (n > 2*m))
   N = (n-1)/m;
   A{1} = [ones(1,N+1)] + [0, repmat(6,[1,N])]*offset;
   for j = 2:m
      A{j} = [j, j-1, j-1] + [0, 7, 1]*offset;
   end
   M(2) = 3*(m-1) + N; % (AA:m-1; A:N+2*(m-1))
   I = [ones([1,m-1]) 0];
else
   N = 1;
   if ((n >= (2*m)-2) & (m > 3))
      N = floor((n+2)/m);
   end
   if N > m % reset offset
      offset = 10^(ceil(log10(N+1))); % label offset
   end
   % Optimize...
   for j = 1:m
      % full cycles, partial forward sequences, partial backward sequences

      % full cycles
      k = 1;
      while (k+m-1) <= length(A{j})
         if A{j}(k) == m
            % replace full cycle
            A{j} = [A{j}(1:k-1), 1+6*offset, A{j}(k+m:end)];
         end
         k = k + 1;
      end

      % multiple full cycles
      f = find(A{j} == 1+6*offset);
      if length(f) > 1 % replace length(f) full cycles
         A{j} = [A{j}(1:f(1)-1), length(f)+6*offset, A{j}(f(end)+1:end)];
      end

      % partial forward sequences
      f = find(A{j} == 1); % find U1 in A{j};
      if f > 1 % substitute Fx for Ux,Ux-1,...,U1
         x = length(find(A{j}(1:f) < offset));
         A{j} = [x+2*offset, A{j}(f+1:end)];
         F(x) = j;
      end

      % partial backwards sequences
      f = find(A{j} == m); % find Um
      if f < length(A{j}) % substitute Bx for Um,Um-1,...,Ux
         x = m+1-length(find(A{j}(f:end) < offset));
         A{j} = [A{j}(1:f-1), x+3*offset];
         B(x) = j;
      end

      % if A{j}(end-1:end) == [Cx B2], replace with [C(x+1) I1]
      if ((floor(A{j}(max(length(A{j}),2)-1)/offset) == 6) & ...
            (A{j}(end) == 2+3*offset))
        
         A{j}(end-1:end) = [A{j}(end-1)+1, 1+offset];
         B(x) = 0;
         I(1) = j;
       end

      % if A{j}(1:2) == [F(m-1) Cx], replace with [Im C(x+1)]
      if ((A{j}(1) == m-1+2*offset) & ...
            (floor(A{j}(min(length(A{j}),2))/offset) == 6))
         A{j}(1:2) = [m+offset, A{j}(2)+1];
         F(x) = 0;
         I(m) = j;
       end

      % near F, near B, J, or K sequences
      if ((length(find(A{j} < offset)) == n) & (n > 2)) % unmodified
         if A{j}(end) == 2 % Fx*I1
            F(A{j}(1)) = j;
            A{j} = [A{j}(1)+2*offset, 1+offset];
            I(1) = j;
         elseif A{j}(1) == m-1 % Im*Bx
            B(A{j}(end)) = j;
            A{j} = [m+offset, A{j}(end)+3*offset];
            I(m) = j;
         elseif (m-A{j}(1)) > (A{j}(end)) % Fz*Jy
            F(A{j}(1)) = j;
            J(A{j}(end)-1) = j;
            A{j} = [A{j}(1)+2*offset, A{j}(end)-1+4*offset];
         else % Kx*Bz
            K(A{j}(1)+1) = j;
            B(A{j}(end)) = j;
            A{j} = [A{j}(1)+1+5*offset, A{j}(end)+3*offset];
         end
      end
   end

   % Remove lowest B and test for improvment vs original.  If no
   % improvment, then use original, and remove next lowest B.
   done = 0;
   f = find(B);
   while (~done & f)
      j = B(f(1));
      B(f(1)) = 0;
      if f(1) == m-1
         A{j} = [A{j}(1:end-1), m, f(1)];
      else
         A{j} = [A{j}(1:end-1), A{j}(end)+1, f(1)];
      end
      if length(A{j}) >= length(original{j})
         for k = 1:length(A{j})
            root = floor(A{j}(k)/offset);
            switch root
            case 1
               I(A{j}(k)-root*offset) = 0;
            case 2
               F(A{j}(k)-root*offset) = 0;
            case 5
               K(A{j}(k)-root*offset) = 0;
            end
         end
         A{j} = original{j};
         done = done - 1;
      end
      done = done + 1;
      if length(f) > 1
         f = f(2:end);
      else
         f = [];
      end
   end

   % Remove highest J and test for improvment vs original.  If no
   % improvment, then use original, and remove next highest J.
   done = 0;
   f = find(J)
   while (~done & f)
      j = J(f(end));
      J(f(end)) = 0;
      if f(end) == 2
         A{j} = [A{j}(1:end-1), 1+offset, 2+offset];
      else
         A{j} = [A{j}(1:end-1), A{j}(end)-1, f(end)+offset];
      end

      if length(A{j}) >= length(original{j})
         for k = 1:length(A{j})
            root = floor(A{j}(k)/offset);
            switch root
            case 1
               I(A{j}(k)-root*offset) = 0;
            case 2
               F(A{j}(k)-root*offset) = 0;
            end
         end
         A{j} = original{j};
         done = done - 1;
      end
      done = done + 1;
      if length(f) > 1
         f = f(1:end-1);
      else
         f = [];
      end
   end

   % Remove lowest K and test for improvment vs original.  If no
   % improvment, then use original, and remove next lowest K.
   done = 0;
   f = find(K);
   while (~done & f)
      j = K(f(1));
      K(f(1)) = 0;
      if f(1) == m-1
         A{j} = [m-1+offset, m+offset, A{j}(2:end)];
      else
         A{j} = [f(1)+offset, A{j}(1)+1, A{j}(2:end)];
      end

      if length(A{j}) >= length(original{j})
         for k = 1:length(A{j})
            root = floor(A{j}(k)/offset);
            switch root
            case 1
               I(A{j}(k)-root*offset) = 0;
            case 3
               B(A{j}(k)-root*offset) = 0;
            end
         end
         A{j} = original{j};
         done = done - 1;
      end
      done = done + 1;
      if length(f) > 1
         f = f(2:end);
      else
         f = [];
      end
   end

   % determine cost (M(2),S) to calculate the F,B,J,K,I, matrices

   % count matrix multiplications

   % check for C(N-1)
   if N > 2
      for j = 1:m
         if any(A{j} == N-1+6*offset)
            C(2) = N-1;
         end
      end
   end

   % check for CN
   if N > 1
      for j = 1:m
         if any(A{j} == N+6*offset)
            C(3) = N;
         end
      end
   end

   % add MMs for CN [and C(N-1)]
   if C(3)
      M(2) = M(2) + N-1; % C(3) (N-1) or C(2) & C(3) [N-2 + 1 = N-1]
   elseif C(2)
      M(2) = M(2) + N-2; % C(2) (N-2)
   end

   % add MMs for A      
   for j = 1:m
      M(2) = M(2) + length(A{j});
   end
   M(2) = M(2) - 1; % correction for length(A{j})-1, C(1) [-m + m-1 = -1]

   % add MMs for B, J, and K sequences
   f = find(B);
   if f
      M(2) = M(2) + m-f(1);
   end
   f = find(J);
   if f
      M(2) = M(2) + f(end)-1;
   end
   f = find(K);
   if f
      M(2) = M(2) + f(end)-1;
   end

   % extra storage; C(1) is not counted since it always has to be stored
   S = length(find([C(2:3), K, J, B, F]));
end
