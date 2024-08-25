% n is the number of intergration points in one dimension
% sdim is the dimension


% w records the weight of each intergration points, is length(w)*1 vector
% x records the local rectangular coordiantes of the intergration points, length(w)*sdim matrix 


function [w,x] = CoorWeight_GL(n,sdim)                           %¸ßË¹º¯Êý
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            % This function determines the abscisas (x) and weights (w)  for the        %
            % Gauss-Legendre quadrature, of order n>1, on the interval [-1, +1].        %
            %   Unlike many publicly available functions, 'GaussLegendre_2' is valid    %
            %   for n>=46. This is due to the fact that 'GaussLegendre_2' does not      %
            %   rely on the build-in Matlab routine 'roots' to determine the roots of   %
            %   the Legendre polynomial, but finds the roots by looking for the         %
            %   eigenvalues of an alternative version of the companion matrix of the    %
            %   n'th degree Legendre polynomial. The companion matrix is constructed    %
            %   as a symmetrical matrix, guaranteeing that all the eigenvalues          %
            %   (roots) will be real. On the contrary, the 'roots' function uses a      %
            %   general form for the companion matrix, which becomes unstable at        %
            %   higher orders n, leading to complex roots.                              %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

            % ?Geert Van Damme
            % geert@vandamme-iliano.be
            % February 21, 2010



            % Building the companion matrix CM
                % CM is such that det(xI-CM)=P_n(x), with P_n the Legendre polynomial
                % under consideration. Moreover, CM will be constructed in such a way
                % that it is symmetrical.
            i   = 1:n-1;
            a   = i./sqrt(4*i.^2-1);
            CM  = diag(a,1) + diag(a,-1);

            % Determining the abscissas (x) and weights (w)
                % - since det(xI-CM)=P_n(x), the abscissas are the roots of the
                %   characteristic polynomial, i.d. the eigenvalues of CM;
                % - the weights can be derived from the corresponding eigenvectors.
            [V,L]   = eig(CM);
            [x1,ind] = sort(diag(L));
            V       = V(:,ind)';
            w1       = 2 * V(:,1).^2;
            
            count = 1;
            if ( sdim == 1 ) 
              x = x1;
              w = w1;

            elseif ( sdim == 2 ) 
              for i = 1:n
                for j = 1:n
                  x(count,:) = [ x1(i), x1(j)];           
                  w(count,1) = w1(i)*w1(j); 
                  count = count+1;
                end
              end

            else % sdim == 3
              for i = 1:n
                for j = 1:n
                  for k = 1:n
                    x(count,:) = [ x1(i), x1(j), x1(k) ];           
                    w(count,1) = w1(i)*w1(j)*w1(k); 
                    count = count+1;
                  end
                end
              end
            end           
            
end