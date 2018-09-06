clear; clc;close all;


delete( "../../1. OUTPUT/test.out")
diary "../../1. OUTPUT/test.out"
#----------------------------------------------
disp("***************************************")
disp("GEOMETRY") 
disp("***************************************")
#----------------------------------------------

disp("---------------------------------------")
disp("Nodal Coordinates: ");
disp("(Node No. - X - Y - Z)");
disp("---------------------------------------")
disp("")
NODES = load('../../0. INPUT/TestModel1.node')
disp("---------------------------------------")
disp("");


disp("---------------------------------------")
disp("Boundary Conditions: 1 = Restrained, 0 = Free")
disp("(Node No. - X? - Y? - Z?)");
disp("---------------------------------------")
disp("")
BOUNDS = load('../../0. INPUT/TestModel1.bound')
disp("---------------------------------------")
disp("");

disp("---------------------------------------")
disp("Element Information:")
disp("(Ele No. - iNode - jNode - Member Type - E - A - I)");
disp("---------------------------------------")
disp("")
ELES = load('../../0. INPUT/TestModel1.ele')
disp("---------------------------------------")
disp("");

disp("---------------------------------------")
disp("Load Information:")
disp("(Node No. - Fx - Fy - Mz)");
disp("---------------------------------------")
disp("")
LOADS = load('../../0. INPUT/TestModel1.load')
disp("---------------------------------------")
disp("");

% specify node coordinates (could only specify non-zero terms)

XYZ = NODES(:,2:3)';
nn = size(XYZ)(2);

% boundary conditions (1 = restrained,  0 = free) (specify only restrained dof's)
BOUN = zeros(nn,3);
for i = 1:size(BOUNDS)(1)
  
  BOUN(BOUNDS(i,1),1:3) = BOUNDS(i,2:4);
  
end

BOUN = BOUN';

% connectivity array
##CON (1,1:2) = [  1   2];   
##CON (2,1:2) = [  2   3];   
##CON (3,1:2) = [  3   4];  
##
##CON = CON'

CON = ELES(:,2:3)';
ne = size(CON)(2);

%Analysis -------------------------------------------------------------------------------

#----------------------------------------------
disp("***************************************")
disp("ANALYSIS") 
disp("***************************************")
disp("")
#----------------------------------------------

for i = 1:ne
    
    beam(i).iNode = CON(1,i);
    beam(i).jNode = CON(2,i);
    
    beam(i).dx = XYZ(1,beam(i).jNode)-XYZ(1,beam(i).iNode);
    beam(i).dy = XYZ(2,beam(i).jNode)-XYZ(2,beam(i).iNode);
    
    beam(i).L = sqrt(beam(i).dx^2+beam(i).dy^2);
    
    beam(i).E = ELES(i,5);
    beam(i).A = ELES(i,6);
    beam(i).I = ELES(i,7);
    
    beam(i).K = getK(beam(i));
    beam(i).K;
    
    beam(i).ae = globalTrans(beam(i).dx/beam(i).L,beam(i).dy/beam(i).L);
    
    beam(i).Ke = (beam(i).ae)'*(beam(i).K)*(beam(i).ae);
    beam(i).Ke;
    
    if i == 1 
        Ka = beam(i).Ke;
    else 
        Ka = blkdiag(Ka,beam(i).Ke);  
    end
end

degreeID = zeros(nn,3);

dof = 1;

for i = 1:nn
    for j = 1:3
        if BOUN(j,i) == 0
           degreeID(i,j) = dof;
           dof = dof + 1;
        else 
           degreeID(i,j) = 0;
        end            
        
    end
end

R = dof;

for i = 1:nn
    for j = 1:3
        if BOUN(j,i) == 1
           degreeID(i,j) = R;
           R = R + 1;
        end          
        
    end
end

for i = 1 : ne
    
    ID(i,1:3) = degreeID(CON(1,i),:);
    ID(i,4:6) = degreeID(CON(2,i),:);
    
end

 dof = dof - 1;
 Af = zeros(ne*6,dof);        %initialize Boolean Colocation matrix
  
 for k = 1:ne               %Populate Boolean colocation matrix
   for id = 1:6
       if ID(k,id) <= dof
            Af( (k-1)*6+id, ID(k,id) ) = 1;
       end
   end
 end
 
 kff = Af'*Ka*Af
 
 % Generate Loads Vector----------------------------------------------------------------
 P = zeros(dof,1);
 nL = size(LOADS)(1);
 
 for i = 1:nL
  
    for j = 1:3
      
      deg = degreeID(LOADS(i,1),j);
      P(deg) = LOADS(i,j+1);
      
    end

 end
 
 %SIMULATE!!!!!!
 
 P
 
 uf = inv(kff)*P
 Qa = Ka*Af*uf;

#----------------------------------------------
disp("***************************************")
disp("Post Process") 
disp("***************************************")
disp("")
#----------------------------------------------

%Nodal Displacements
NodalDisp = zeros(nn,4);

NodalDisp(:,1) = NODES(:,1);

for i = 1:size(NODES)(1)

  for j = 1:3
    
    if degreeID(i,j) > dof 
      
    else
      NodalDisp(i,j+1) = uf(degreeID(i,j));
    end
    
    
  end
 
end


disp("---------------------------------------")
disp("Nodal Displacements: ");
disp("(Node No. - Tx - Ty - Rz)");
disp("---------------------------------------")
disp("")
NodalDisp
disp("---------------------------------------")

%Member end forces (Global):
GlobalForces = zeros(ne,7);
GlobalForces(:,1) = ELES(:,1);
for i = 1:ne
  
  GlobalForces(i,2:7) = Qa(1 + (i-1)*6 : 6 + (i-1)*6);

end

disp("---------------------------------------")
disp("Global Member Forces: ");
disp("(Element No. - iFx - iFy - iMz - jFx - jFy - jMz)");
disp("---------------------------------------")
disp("")
GlobalForces
disp("---------------------------------------")

%Member end forces (GLOBAL):

LocalForces = zeros(ne,7);
LocalForces(:,1) = ELES(:,1);

for i = 1:ne
  
  LocalForces(i,2:7) = beam(i).ae*GlobalForces(i,2:7)';
  
end

disp("---------------------------------------")
disp("Local Member Forces: ");
disp("(Element No. - iFx - iFy - iMz - jFx - jFy - jMz)");
disp("---------------------------------------")
disp("")
LocalForces
disp("---------------------------------------")

%Support Reactions

Support = zeros(nn,4);

Support(:,1) = NODES(:,1);

%for each node
for node = 1:nn
  
    %for each degree of freedom
    for degree = 1:3
      
        if degreeID(node,degree) > dof
          
          R = 0;
          for element = 1:ne
            
            iNode = CON(1,element);
            jNode = CON(2,element);
            
            if iNode == node
              
              R = R + GlobalForces(element,degree+1);
              
            end
            
            if jNode == node
              
              R = R + GlobalForces(element,degree+4);
              
            end
            
          end %for each element
          
          Support(node,degree+1) = R;
              
        end %if degree is restrained
        
    end %for each degree
 
end %for each node
disp("---------------------------------------")
disp("Nodal Support Reactions: ");
disp("(Node No. - Fx - Fy - Mz)");
disp("---------------------------------------")
Support
sumFx = sum(Support(:,2))
sumFy = sum(Support(:,3))
sumMz = sum(Support(:,4))
disp("---------------------------------------")


diary off;
