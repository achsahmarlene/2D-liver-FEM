%{

We are taking cross-section of a tumor which is circular in shape. 
133 nodes are used to model the circular tumor. 
We are taking a series of dodecagons(12 sided polygon) where the vertices will
act as the nodes.

%}

% finding the node points

r1 = 3/2:-0.5/2:0.5/2;
r1 = reshape(r1, size(r1,2), 1);
theta1 = 0:pi/6:(2*pi - pi/6);
xp1 = r1*cos(theta1);
xp1 = xp1';
xp1 = xp1(:);
yp1 = r1*sin(theta1);
yp1 = yp1';
yp1 = yp1(:);

r2 = 2.75/2:-0.5/2:0.75/2;
r2 = reshape(r2, size(r2,2), 1);
theta2 = pi/12:pi/6:(2*pi - pi/12);
xp2 = r2*cos(theta2);
xp2 = xp2';
xp2 = xp2(:);
yp2 = r2*sin(theta2);
yp2 = yp2';
yp2 = yp2(:);

nodes = [xp1 yp1; xp2 yp2; 0 0];
nodes = round(nodes, 6); 
% the node matrix is of size (133,2) where the 1st column corresponds to
% x-coordinate and 2nd column corresponds to y-coordinate of a node 

cm = zeros(252, 3); % connectivity matrix

for k=0:4
    for i=1:12
        m1 = mod(i,12);
        m2 = mod(i+1,12);
        if m1==0
            m1=12;
        end
        if m2==0
            m2=12;
        end
        m1 = m1 + k*12;
        m2 = m2 + k*12;
        cm(4*(i+k*12) -3, :) = [m1 m2 i+k*12+72];
        cm(4*(i+k*12) -2, :) = [m2 m2+12 i+k*12+72];
        cm(4*(i+k*12) -1, :) = [m2+12 m1+12 i+k*12+72];
        cm(4*(i+k*12), :) = [m1+12 m1 i+k*12+72];
    end
end

for i=1:12
    m1 = mod(i,12);
    if m1==0
        m1=12;
    end
    m2 = mod(i+1,12);
    if m2==0
        m2=12;
    end
    cm(240 + i, :) = [60+m1 60+m2 133];
end
%{
we found the connectivity matrix of all the element
%}

K_global = zeros(133,133); % global stiffness matrix
k = 0.57; % thermal conductivity of liver tumor tissue

for i=1:252
    K_local = zeros(3,3); % local stiffness matrix
    N = ones(3,3); % shape matrix
   
    N(1,2) = nodes(cm(i,1), 1);
    N(1,3) = nodes(cm(i,1), 2);
    N(2,2) = nodes(cm(i,2), 1);
    N(2,3) = nodes(cm(i,2), 2);
    N(3,2) = nodes(cm(i,3), 1);
    N(3,3) = nodes(cm(i,3), 2);

    area = (1/2)*det(N); % area of an element
    B = zeros(2,3);

    B(1,1) = N(2,3) - N(3,3);
    B(2,1) = N(3,2) - N(2,2);
    B(1,2) = N(3,3) - N(1,3);
    B(2,2) = N(1,2) - N(3,2);
    B(1,3) = N(1,3) - N(2,3);
    B(2,3) = N(2,2) - N(1,2);
    B = (1/(2*area)).*B; 

    K_local(:,:) = B'*B;
    K_local(:,:) = (k*area).*K_local(:,:); % local stiffness matrix

    % adding the values from local stiffness matrix to global stiffness
    % matrix at their corresponding positions
    K_global(cm(i,1),cm(i,1)) = K_global(cm(i,1),cm(i,1)) + K_local(1,1); 
    K_global(cm(i,1),cm(i,2)) = K_global(cm(i,1),cm(i,2)) + K_local(1,2);
    K_global(cm(i,1),cm(i,3)) = K_global(cm(i,1),cm(i,3)) + K_local(1,3);
    K_global(cm(i,2),cm(i,1)) = K_global(cm(i,2),cm(i,1)) + K_local(2,1);
    K_global(cm(i,2),cm(i,2)) = K_global(cm(i,2),cm(i,2)) + K_local(2,2);
    K_global(cm(i,2),cm(i,3)) = K_global(cm(i,2),cm(i,3)) + K_local(2,3);
    K_global(cm(i,3),cm(i,1)) = K_global(cm(i,3),cm(i,1)) + K_local(3,1);
    K_global(cm(i,3),cm(i,2)) = K_global(cm(i,3),cm(i,2)) + K_local(3,2);
    K_global(cm(i,3),cm(i,3)) = K_global(cm(i,3),cm(i,3)) + K_local(3,3);
end

K_global = round(K_global,4);

% heat supplied by trocar at all elements
% only the center is supplied with heat by trocar
trocar_heat = zeros(133,1);
trocar_heat(133,:) = 10;

% metabolic heat is uniform across the cross-section         
% the metabolic heat at each element turns out to be 0.0018
metabolic_heat = zeros(133,1) + 0.0018;

% inverse of global stiffness matrix
K_global_inv = inv(K_global);

% due to symmetry, we can say that all boundary points have same reaction
% this assumption reduces our problem from 12 linear equations to just one
% linear equation to solve
R = (50 - K_global_inv(1,:)*(trocar_heat + metabolic_heat))/sum(K_global_inv(1,1:12),'all');
% we got reaction as -0.0199

fprintf('Value of reaction is %f\n',R);

% reaction at each node
% reaction at all nodes except boundary nodes is zero
reaction = zeros(133,1);
reaction(1:12,:) = zeros(12,1) + R;

% finding the total heat
heat = trocar_heat + metabolic_heat + reaction;
heat = reshape(heat, size(heat,1), 1);

% temperature matrix to find temperature at each node
temp = zeros(133,1);
for i=1:133
    temp(i,:) = K_global_inv(i,:)*heat;
end

% displaying temperature required to be supplied at the center of the tumor
fprintf('Value of temperature at the center of the tumor is %f\n',temp(133,:));