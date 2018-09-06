function [ K ] = getK( beam)
  E = beam.E;
  A = beam.A;
  I = beam.I;
  L = beam.L;

  K = [A*L^2/I        0        0       -A*L^2/I      0       0;
       0              12       6*L        0         -12      6*L;
       0             6*L     4*L^2        0        -6*L     2*L^2;
       -A*L^2/I       0         0       A*L^2/I       0       0;
       0             -12       -6*L       0         12      -6*L;
       0              6*L      2*L^2      0        -6*L      4*L^2];
       
       
  K = E*I/L^3 * K;

end

