%% Written by Dipjyoti Das, MCDB, Yale
clear all; close all;
 
 
%% Read the file.
filename1 = '3D_cell_motion.dat'; 

zmid=5*0.5;

his1 = dlmread(filename1); 



 T = his1(:,1); % Time Step
 X = his1(:,2); % Second column is the x position of the particles
 Y = his1(:,3); % 3rd column is the y position of the particles
 Z = his1(:,4); % 4th column is the z position of the particles
 
 VX = his1(:,5); % 5th column is the x component of  particle velocity
 VY = his1(:,6); % 6th column is the y component of  particle velocity
 VZ = his1(:,7); % 7th column is the z component of  particle velocity
 
 Index = his1(:,8); % 8th column is the particle index
 
 TgState = his1(:,9); % 9th column is the transgenic state of particles
 yPZ = his1(:,10);
 
 
 
ymax=max(Y);
ymin=min(Y);
zmax=max(Z);
zmin=min(Z);
xmax=max(X);
xmin=min(X);

tmin = min(T); % Time at beginning.
tmax = max(T); % Time at end
dt = 1;  % How often would you like to capture the frame.
 %% Define movie name etc.
 
 figure(1), % Shall open an empty figure window. Adjust it to the size you want. 
 % Ideally, you should have a fixed figure with specified position and size, but 
 % I find that very constrained. So, I just open an empty figure window, and adjust 
 % it according to my need. Don't touch it until the code completes. It'll interrupt 
 % the program, and the movie won't be produced. I sometime use this method to kill a 
 % process and redo the movie to my like.   
 
% clear movObj
% movObj = VideoWriter('Particle_Evolution.avi'); % Give whatever name you want to give. Can be automated. 
%movObj.FrameRate = 10; % Change to whatever value you want. 
 %.Quality = 75; % Any number between 0-100. 
 
 %(movObj);
 
 %% The best part: Frame capture 
 
 for t = tmin:dt:tmax 
	id = T == t; % All data index for time t. 
	 yy=yPZ(id);
     yUP=yy(1)+2.0 ;
    
    id2= TgState==1 & T==t;
    
    id3= T==t & VY<0 & Y<=yUP & Y>4 & Z>=zmid;
    
    
	x = X(id); 
	y = Y(id);
    z = Z(id);
	vx = VX(id); 
	vy = VY(id); 
	vz = VZ(id);
    
    xp = X(id2); 
	yp = Y(id2);
    zp = Z(id2);
 
    xpp = X(id3); 
	ypp = Y(id3);
    zpp = Z(id3);
    
  plot3(x,y,z,'ob','LineWidth',1.2,'MarkerSize',12,'MarkerFaceColor','c')
  hold on;
  plot3(xp,yp,zp,'ob','LineWidth',1.2,'MarkerSize',12,'MarkerFaceColor','g')
  plot3(xpp,ypp,zpp,'ob','LineWidth',1.2,'MarkerSize',12,'MarkerFaceColor','r')
   %quiver3(x,y,z,vx,vy,vz,0,'k','LineWidth',1.4) % Plot the velocity vectors
	%axis square  % add other axis properties if you feel like. 
	xlabel('X'); 
	ylabel('Y');
    zlabel ('Z')
	title(['Time = ' num2str(t)])
    
    axis equal
     
  
 
 xlim([xmin-1 xmax+1]);
 ylim([ymin ymax+1]);
 zlim([zmin zmax+0.1]); 
 % view(90,  90 ); %%%% Use this for dorsal (top) view of the tailbud
  view(90,  0 );   %%%% Use this for lateral (side) view of the tailbud
 
 hold off;
 
 
 %axis vis3d
 
	% Movie stuff	
    M(t)=getframe;
	frame = getframe(gcf); % Capture the frame
	%writeVideo(movObj,frame); % Add the frame to the movie.
 end
 
 
 numtimes=1;
 fps=29;
 movie(M,numtimes,fps)
 %close(movObj); % Close the movie variable and write the file. 
%movie2avi(M, 'particle-history.avi'); %, 'compression', 'none');





 
 
 
 %plot3(XPSMR,YPSMR,ZPSMR,'ob','LineWidth',1.2,'MarkerSize',12,'MarkerFaceColor','b');
 
 
 
