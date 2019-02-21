% some svd figure things out stuff

At =[1, 3; 3 10];

[u,s,v] = svd(At)


I = eye(2);

vi = v'*I;
svi = s*vi;
usvi = u*svi;
AI = At*I
% lets look at some some plots
figure(1), 
plot([0 I(1,1)], [0 I(2,1)], 'b'), hold on
plot([0 I(1,2)], [0 I(2,2)], 'b:')
grid on

plot([0 vi(1,1)], [0 vi(2,1)], 'r'), hold on
plot([0 vi(1,2)], [0 vi(2,2)], 'r:'), 

plot([0 svi(1,1)], [0 svi(2,1)], 'g'), hold on
plot([0 svi(1,2)], [0 svi(2,2)], 'g:'), 
 
plot([0 usvi(1,1)], [0 usvi(2,1)], 'y'), hold on
plot([0 usvi(1,2)], [0 usvi(2,2)], 'y:'), 
   
plot([0 AI(1,1)], [0 AI(2,1)], 'k--'), hold on
plot([0 AI(1,2)], [0 AI(2,2)], 'k:'), 


%% First singular values and a vector
x = [1; 0];

vx = v'*x;
svx = s*vx;
usvx = u*svx;
Ax = At*x;
figure(2)
plot([0 x(1)],[0 x(2)],'b'), hold on, grid on
plot([0 vx(1)],[0 vx(2)], 'r')
plot([0 svx(1)],[0 svx(2)],'g')
plot([0 usvx(1)],[0 usvx(2)],'y')
plot([0 Ax(1)],[0 Ax(2)], 'k--')