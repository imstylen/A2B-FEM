function [ K ] = getK( beam )
  typ = beam.typ; 
  E = beam.E;
  A = beam.A;
  I = beam.I;
  L = beam.L;

if typ == 0

  K = [A*L^2/I        0        0       -A*L^2/I      0       0;
       0              12       6*L        0         -12      6*L;
       0             6*L     4*L^2        0        -6*L     2*L^2;
       -A*L^2/I       0         0       A*L^2/I       0       0;
       0             -12       -6*L       0         12      -6*L;
       0              6*L      2*L^2      0        -6*L      4*L^2];
       
       
  K = E*I/L^3 * K;
  
  
elseif typ == 1
  K = [A*L^2/I        0        0       -A*L^2/I      0         0;
       0              3        0         0          -3        3*L;
       0              0        0         0           0         0;
       -A*L^2/I       0        0       A*L^2/I       0         0;
       0             -3        0         0           3        -3*L;
       0              3*L      0         0          -3*L      3*L^2];
       
       
  K = E*I/L^3 * K;
  
elseif typ == 2
  
  K = [A*L^2/I        0         0       -A*L^2/I      0      0;
       0              3         3*L         0        -3      0;
       0              3*L      3*L^2        0       -3*L     0;
       -A*L^2/I       0         0       A*L^2/I       0      0;
       0             -3        -3*L         0         3      0;
       0              0         0           0         0      0];
       
       
  K = E*I/L^3 * K;
  
elseif typ == 3
  
  K = [1        0        0       -1      0       0;
       0        0        0        0      0       0;
       0        0        0        0      0       0;
       -1       0        0        1      0       0;
       0        0        0        0      0       0;
       0        0        0        0      0       0];
       
       
  K = E*A/L * K; 

end


  

end

