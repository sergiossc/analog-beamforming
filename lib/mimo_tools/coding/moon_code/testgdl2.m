% Test the generalized distributive law, as implemented on 
% Junction Graphs

% Copyright 2004 by Todd K. Moon
% Permission is granted to use this program/data
% for educational/research only


% Describe the Tanner graph for the code
%       1 2 3 4 5 6 7 8 9 10
% 11   [1 1 1 0 0 1 1 0 0 1
% 12    1 0 1 0 1 1 0 1 1 0
% 13    0 0 1 1 1 0 1 0 1 1
% 14    0 1 0 1 1 1 0 1 0 1
% 15    1 1 0 1 0 0 1 1 1 0];

% Desriptions of interconnections
C{1} = [1 11 12 15];   % variable to check nodes
C{2} = [2 11 14 15];
C{3} = [3 11 12 13];
C{4} = [4 13 14 15];
C{5} = [5 12 13 14];
C{6} = [6 11 12 14];
C{7} = [7 11 13 15];
C{8} = [8 12 14 15];
C{9} = [9 12 13 15];
C{10} = [10 11 13 14];
C{11} = [11 1 2 3 6 7 10];   % check nodes to variable nodes
C{12} = [12 1 3 5 6 8 9];
C{13} = [13 3 4 5 7 9 10];
C{14} = [14 2 4 5 6 8 10];
C{15} = [15 1 2 4 7 8 9];

% list of local kernal functions
dens = 'fyx0';  % description of the memoryless channel model

F{1} = dens;
F{2} = dens;
F{3} = dens;
F{4} = dens;
F{5} = dens;
F{6} = dens;
F{7} = dens;
F{8} = dens;
F{9} = dens;
F{10} = dens;
F{11} = 'parity';
F{12} = 'parity';
F{13} = 'parity';
F{14} = 'parity';
F{15} = 'parity';

% List of variables at each node
S{1} = 1;
S{2} = 2;
S{3} = 3;
S{4} = 4;
S{5} = 5;
S{6} = 6;
S{7} = 7;
S{8} = 8;
S{9} = 9;
S{10} = 10;
S{11} = [1 2 3 6 7 10];
S{12} = [1 3 5 6 8 9];
S{13} = [3 4 5 7 9 10];
S{14} = [2 4 5 6 8 10];
S{15} = [1 2 4 7 8 9];

% Received data
% A used in chapter
A = [
1 1 1 0 0 1 1 0 0 1
1 0 1 0 1 1 0 1 1 0
0 0 1 1 1 0 1 0 1 1
0 1 0 1 1 1 0 1 0 1
1 1 0 1 0 0 1 1 1 0];

Apinv = [1 0 1 1 0;
	     0 1 1 0 1;
		 0 1 0 1 1;
		 1 1 0 1 0;
		 1 0 1 0 1];
H = mod(Apinv*A,2);   % [I P]

[M,N] = size(A);
K = N-M;
P = H(:,N-K+1:N);
G = [P; eye(K)];  % now A*G = 0 (mod 2)
m = [1 0 1 0 1]';
c = mod(G*m,2);
t = 2*(2*c-1);
p1 = [.22  .16  .19  .48 .55  .87 .18 .79 .25 .76]; 
a = 2; sigma2 = 2;
y =  log((1./p1)-1)/(-2*a)*sigma2;

Y{1} = [y(1) sigma2 a];   % observation and variance
Y{2} = [y(2) sigma2 a];   % observation and variance
Y{3} = [y(3) sigma2 a];   % observation and variance
Y{4} = [y(4) sigma2 a];   % observation and variance
Y{5} = [y(5) sigma2 a];   % observation and variance
Y{6} = [y(6) sigma2 a];   % observation and variance
Y{7} = [y(7) sigma2 a];   % observation and variance
Y{8} = [y(8) sigma2 a];   % observation and variance
Y{9} = [y(9) sigma2 a];   % observation and variance
Y{10} = [y(10) sigma2 a];   % observation and variance
Y{11} = [];
Y{12} = [];
Y{13} = [];
Y{14} = [];
Y{15} = [];

Niter = 10;
Nin = 10;   % number of nodes to compute belief
Q = 2;  % binary values
Nvar = 10;
gdl