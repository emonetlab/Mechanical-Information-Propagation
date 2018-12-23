!!!!!!!!!! Written by Dipjyoti Das, Holley Lab, MCDB, Yale, 2018
!!!!! see Das, Jülich and Schwendinger-Schreck et al., 2018, "Organization of embryonic morphogenesis via mechanical information"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Local variables and parameters !!!!!!!!!!!!!!!!!!!!!!!

MODULE numz
  IMPLICIT NONE
  INTEGER,PARAMETER:: DP=KIND(1.0d0)
  REAL(DP),PARAMETER:: pi=3.14159265358979_DP
END MODULE numz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE constants
 USE numz
 IMPLICIT NONE
REAL(DP),PARAMETER:: L=10.0_DP  !! initial length
REAL(DP),PARAMETER:: R=5.0_DP !L*0.5 !!/2.0_DP  !! Radius of PZ (=L/2)
REAL(DP),PARAMETER:: mu=1.0_DP  !! mobility
REAL(DP),PARAMETER:: Vo=1.0_DP  !! magnitude of self-propulsion velocity
REAL(DP),PARAMETER:: tau=1.0_DP  !! relaxation time for the alingment of self-propulsion vectors

REAL(DP),PARAMETER:: Fadh=1.000_DP !! Maximums of the Replusive and adhesive forces
REAL(DP),PARAMETER:: ao=(5.0_DP/6.0_DP)*0.5_DP !! average radius of cells
REAL(DP),PARAMETER:: pack_frac=0.95_DP  !! packing fraction of cells (initially)
REAL(DP),PARAMETER:: Vcell=(4.0_DP/3.0_DP)*pi*ao**3 !!average cell volume

!!!!%%%%%% initial numbers of cells in ADM/NT,DM, TO/PZ, PSM (all integers) %%%%%%%%%%%%
!!!! This determines how many cells are inside each region (initially)

REAL(DP),PARAMETER::zmid= R*0.500_DP !!Hight of the mid-line from z=0 dividing NT and PSM
INTEGER,PARAMETER:: NPtop=  NINT( (1.0/4.0)*((4.0_DP/3.0_DP)*pi*R**3)*(pack_frac/Vcell) ) 
INTEGER,PARAMETER:: NPneck=  NINT( (0.5*pi*R**2)*(L/5.0)*(pack_frac/Vcell) ) 
INTEGER,PARAMETER:: NADM= NINT((R**2/2.0)*( (2.0*acos(zmid/R)) - sin(2.0*acos(zmid/R)))*L *(pack_frac/Vcell)) !!!!!!!!! area of circular sector= (R**2/2)(theta - sin(theta))
INTEGER,PARAMETER:: NPSM=   NINT( ((0.5*pi*R**2)- (R**2/2.0)*( (2.0*acos(zmid/R)) - sin(2.0*acos(zmid/R)) ) )*L*(pack_frac/Vcell) )  
INTEGER,PARAMETER:: Ntot_initial = NPtop+NADM+NPSM+NPneck
INTEGER,PARAMETER:: Narray=Ntot_initial*1000

!!!!! time & sampling %%%%%%%%%%%%%%%%%%%%%%%%
INTEGER,PARAMETER:: Maxstep=9000   !!total simulation time (in integral steps)
REAL(DP),PARAMETER:: h=0.005     !!time-increment
INTEGER,PARAMETER:: sample_step=30  !! after how many step-intervals, we are goin g to sample the system
INTEGER,PARAMETER:: Ttrans=3000 !100*sample_step

!!!!! CELL INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%

INTEGER,PARAMETER:: step_input = 5 !! Step #, after which 1 cell is introduced (denoted as gamma in our paper)
END MODULE constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!This fortran 90 script simulate the 3D cell motion following a Vicsek dynamics in a zebrafish tailbud
!!The tailbud geometry is a half-cylinder with rigid boundaries
!! For details see D. Das et al., "Organization of embryonic morphogenesis via mechanical information", 2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM horseshoe
 USE numz
 USE constants 

 IMPLICIT NONE 

!!!!!!! Difine variables !!!!!!!!!!!

!!!!%%%%%%%%%%%%%%%%%% position and velocity arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(DP),DIMENSION(1:Narray, 1:3)::position, velocity !!position (x,y,z) & velocity (Vx, Vy, Vz) of each particle
REAL(DP),DIMENSION(1:Narray,1:3)::self_prop !! self-propulsion direction of each particle
REAL(DP),DIMENSION(1:Narray)::radi !! radius of each particle
REAL(DP),DIMENSION(1:Narray, 1:3):: position_prev


REAL(DP),DIMENSION(1:Narray):: Frep !!! Repulsion on each particle
INTEGER,DIMENSION(1:Narray):: Tg_state !!! 0=normal, 1=transgenic

REAL(DP)::eta !!!!!!!! scalr noise strength running between [0,1]
!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!! sampling %%%%%%%%%%%%%%%%%%%%%%%%
INTEGER:: sample_counter=1
REAL(DP),DIMENSION(1:3):: e_cap, e_cap_temp,vel_cap
REAL(DP),DIMENSION(1:3,1:3)::Rot_matrix,Identity
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!%%%%% variables for the tailbud geometry %%%%%%%%
REAL(DP):: ymax,yPZ,ymin   

!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!

!!!%%%%%%%% LOCAL variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INTEGER::idum,k,Ntot,step,i,j,n,Ntemp,dN, nn
REAL(DP):: ran2,rc,polar_angle,azimuth_angle,rad, testx, testz, anglemin, anglemax
REAL(DP):: volume, volume_temp,dvol,dl
REAL(DP):: fxsum,fysum,fzsum,fx,fy,fz,fwallx,fwally,fwallz
REAL(DP):: xtemp,ytemp,ztemp,Rot_angle, noise_angle
REAL(DP):: vx0,vy0,vz0,ux0,uy0,uz0
REAL(DP):: xx1,yy1,zz1,xx2,yy2,zz2 
REAL(DP):: Neighbor, dnn, yADMtop,yADM_min


!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!! Initialize all arrays 
fwallx=0.0_DP
fwally=0.0_DP
fwallz=0.0_DP

position=0.0_DP
velocity=0.0_DP
radi=0.0_DP
self_prop=0.0_DP
position_prev=0.0_DP
Identity = 0.0_DP
   DO k=1,3
     Identity(k,k) = 1.000_DP
   END DO 
e_cap = 0.0_DP         
e_cap_temp = 0.0_DP   
vel_cap = 0.0_DP
Rot_matrix = 0.0_DP            




eta=0.700_DP !!! Initialize noise 
Frep=30.0_DP !!! Initialize repulsion
Tg_state=0	 !!! Initialize transgenic state


!!!!! Build the initial horse-shoe geometey!!!!!!!!!!!!!!
!ymax=L+R   !+L/30.0_DP
ymax=L+(L/5.0)
yPZ=L
ymin=0.0_DP
yADM_min=L/2.0

!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPEN(unit=21,file='3D_cell_motion.dat',status='unknown')


CALL SYSTEM_CLOCK(COUNT=idum)
!!!!!****** INITIALLY CELLS TAKEN RANDOMLY IN DIFFERENT REGIONS*************

!!!! ADM region : put NADM number of cells
     do k=1,NADM
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   
	    anglemin = (pi/2.0) - (acos(zmid/R))
	    anglemax = ( (pi/2.0) - (acos(zmid/R)) ) + 2.0*(acos(zmid/R))
	    
    50 polar_angle= anglemin +(anglemax-anglemin)*ran2(idum)   !! 30 to 150 degree angle
	   rad=(R-rc)*dsqrt(ran2(idum))
	   
	   testx= rad*dcos(polar_angle)
	   testz= rad*dsin(polar_angle)
	  
	 IF (testz >= (zmid+rc) ) then 
	   position(k,1)= testx
	   position(k,3)= testz
	   position(k,2)= ymin+rc + (yPZ-ymin-2.0_DP*rc)*ran2(idum)   
	 ELSE
	  GOTO 50
	 END IF
	 
	 azimuth_angle=2.0_DP*pi*ran2(idum)
	 polar_angle= dacos(2.0*ran2(idum) - 1.0)	!!!!!!!!!! cos^-1[-1,1]   
	 self_prop(k,1)= dsin(polar_angle)*dcos(azimuth_angle)
	 self_prop(k,2)= dsin(polar_angle)*dsin(azimuth_angle)  
	 self_prop(k,3) = dcos(polar_angle)
	 
	   position_prev(k,1)= position(k,1)
	   position_prev(k,2)= position(k,2)
	   position_prev(k,3)= position(k,3)
	end do
	

	
!!! PZ+DM Region: put NPtop number of cells

	do k=(NADM+1),(NADM+NPtop)
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	51 azimuth_angle=pi*ran2(idum)
	   polar_angle= dacos(ran2(idum)) !!!!!!!!!!!!! cos^(-1)[0,1]
	   rad= (R-rc)*(ran2(idum))**(1.0_DP/3.0_DP)
	   
	   testx=rad*dsin(polar_angle)*dcos(azimuth_angle)
	   testz=rad*dsin(polar_angle)*dsin(azimuth_angle)
	   
	  IF ( testz>=rc ) then 
	   position(k,1)= testx
	   position(k,3)= testz
	   position(k,2)= ymax+rc + rad*dcos(polar_angle) 
	   ELSE
	  GOTO 51
	 END IF  
	 
	 azimuth_angle=2.0_DP*pi*ran2(idum)
	 polar_angle= dacos(2.0*ran2(idum) - 1.0)	   
	 self_prop(k,1)= dsin(polar_angle)*dcos(azimuth_angle)
	 self_prop(k,2)= dsin(polar_angle)*dsin(azimuth_angle)  
	 self_prop(k,3) = dcos(polar_angle) 
	   
	   position_prev(k,1)= position(k,1)
	   position_prev(k,2)= position(k,2)
	   position_prev(k,3)= position(k,3)
	end do	
	
	
!!! PSM-left Region: put NPSM/2 number of cells

   do k=(NADM+NPtop+1),NINT(NADM+NPtop+NPSM/2.0)
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   
   52  polar_angle= (pi/2.0) +(pi-(pi/2.0))*ran2(idum)
	   rad=(R-rc)*dsqrt(ran2(idum))
	   testx= rad*dcos(polar_angle)
	   testz= rad*dsin(polar_angle)
	  
	 IF (testz <= (zmid-rc) .and. testz>=rc .and. testx <= -rc) then 
	   position(k,1)= testx
	   position(k,3)= testz
	   position(k,2)= ymin+rc + (yPZ-ymin-2.0_DP*rc)*ran2(idum) 
	   
	 ELSE
	  GOTO 52
	 END IF
	   
	 azimuth_angle=2.0_DP*pi*ran2(idum)
	 polar_angle= dacos(2.0*ran2(idum) - 1.0)	   
	 self_prop(k,1)= dsin(polar_angle)*dcos(azimuth_angle)
	 self_prop(k,2)= dsin(polar_angle)*dsin(azimuth_angle)  
	 self_prop(k,3) = dcos(polar_angle) 
	  
	   position_prev(k,1)= position(k,1)
	   position_prev(k,2)= position(k,2)
	   position_prev(k,3)= position(k,3)
	end do
	
!!! PSM-right Region: put NPSM/2 number of cells

   do k=1+NINT(NADM+NPtop+NPSM/2.0),(NADM+NPtop+NPSM)
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	53  polar_angle=(pi/2.0)*ran2(idum)
	   rad=(R-rc)*dsqrt(ran2(idum))
	   testx= rad*dcos(polar_angle)
	   testz= rad*dsin(polar_angle)
	  
	 IF (testz <= (zmid-rc) .and. testz>=rc .and. testx >= rc) then 
	   position(k,1)= testx
	   position(k,3)= testz
	   position(k,2)= ymin+rc + (yPZ-ymin-2.0_DP*rc)*ran2(idum) 
	   
	 ELSE
	  GOTO 53
	 END IF
	   
	 azimuth_angle=2.0_DP*pi*ran2(idum)
	 polar_angle= dacos(2.0*ran2(idum) - 1.0)	   
	 self_prop(k,1)= dsin(polar_angle)*dcos(azimuth_angle)
	 self_prop(k,2)= dsin(polar_angle)*dsin(azimuth_angle)  
	 self_prop(k,3) = dcos(polar_angle) 
	 
	   position_prev(k,1)= position(k,1)
	   position_prev(k,2)= position(k,2)
	   position_prev(k,3)= position(k,3)
	end do

!!!! ADM-PZ-neck region : put NPneck number of cells

     do k=(NADM+NPtop+NPSM+1),Ntot_initial
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   
	54  polar_angle= pi*ran2(idum)
	   rad=(R-rc)*dsqrt(ran2(idum))
	   testx= rad*dcos(polar_angle)
	   testz= rad*dsin(polar_angle)
	  
	 IF (testz >= rc ) then 
	   position(k,1)= testx
	   position(k,3)= testz
	   position(k,2)= yPZ+rc + (ymax-yPZ-rc)*ran2(idum) 
	   
	 ELSE
	  GOTO 54
	 END IF
	 
	 azimuth_angle=2.0_DP*pi*ran2(idum)
	 polar_angle= dacos(2.0*ran2(idum) - 1.0)	   
	 self_prop(k,1)= dsin(polar_angle)*dcos(azimuth_angle)
	 self_prop(k,2)= dsin(polar_angle)*dsin(azimuth_angle)  
	 self_prop(k,3) = dcos(polar_angle) 
	 
	   position_prev(k,1)= position(k,1)
	   position_prev(k,2)= position(k,2)
	   position_prev(k,3)= position(k,3)
	end do 
!!!!!***********************************************************************

!do k=1,Ntot_initial
!write(21,*) position(k,1),position(k,2),position(k,3),self_prop(k,1),self_prop(k,2),self_prop(k,3),k
!end do


!!! INITIALIZATION !!!!!!!!!!!!!

volume= (0.5_DP*pi*R**2)*ymax + 0.25_DP*(4.0_DP/3.0_DP)*pi*R**3
 Ntot=Ntot_initial
sample_counter=1



DO step=1,Maxstep !!!! time-loop starts



   Do i=1,Ntot   !!!! BEGIN Cell Index Loop
   
   
	
        CALL SYSTEM_CLOCK(COUNT=idum)

         xtemp=position(i,1)
         ytemp=position(i,2)
  		 ztemp=position(i,3)
         rc=radi(i)
         
    !!!*********** INTRODUCE TG cells at PZ/TO after Ttrans******************
    
    if( step>=Ttrans .and. position(i,2) >= yPZ .and. position(i,3)<=zmid )then
    if(ran2(idum)<=0.8 .and. Tg_state(i)== 0)then   
     Tg_state(i)=1
     end if
    end if
    !!****************************************************************   

   !!Calculate the force on ith cell due to all others **********************
       fxsum=0.0_DP
       fysum=0.0_DP
	   fzsum=0.0_DP
	   
         do j=1,Ntot   !!! BEGIN: Sum the forces Loop
         
         xx1=position(i,1)
         yy1=position(i,2)
         zz1=position(i,3)
         xx2=position(j,1)
         yy2=position(j,2)
         zz2=position(j,3)	         

call cell_force(yPZ,xx1,yy1,zz1,xx2,yy2,zz2,radi(i),radi(j),Tg_state(j),fx,fy,fz)	
	         fxsum=fxsum+fx
	         fysum=fysum+fy
	         fzsum=fzsum+fz
	
        end do        !!!!END: Sum the forces Loop
        
     
        	 
     !!!!! UPDATE positions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      position(i,1) = position(i,1) + (Vo*self_prop(i,1)+ mu*fxsum)*h
      position(i,2) = position(i,2) + (Vo*self_prop(i,2)+ mu*fysum)*h
      position(i,3) = position(i,3) + (Vo*self_prop(i,3)+ mu*fzsum)*h
      
 
  CALL boundary(ymax,yPZ,rc,xtemp,ytemp,ztemp,position(i,1),position(i,2),position(i,3)) !!! reflecting boundary     
       
    

     !!	******************UPDATE THE SELF-PROPULSION DIRECTION************************
	      velocity(i,1)=( position(i,1)-xtemp)/h
	      velocity(i,2)=( position(i,2)-ytemp)/h
	      velocity(i,3)=( position(i,3)-ztemp)/h
	  vel_cap(1)=velocity(i,1)/(dsqrt(velocity(i,1)**2+velocity(i,2)**2+velocity(i,3)**2))
	  vel_cap(2)=velocity(i,2)/(dsqrt(velocity(i,1)**2+velocity(i,2)**2+velocity(i,3)**2))
	  vel_cap(3)=velocity(i,3)/(dsqrt(velocity(i,1)**2+velocity(i,2)**2+velocity(i,3)**2))
	  
	  e_cap_temp=(/0.00_DP, vel_cap(3), -vel_cap(2)/)  !!!!!!! (v^ X x^)
	  e_cap_temp = e_cap_temp/NORM2(e_cap_temp)
	  
	 
	  
	  Rot_matrix(1,:)=(/ 0.0_DP, -vel_cap(3), vel_cap(2)/)
	  Rot_matrix(2,:)=(/ vel_cap(3), 0.0_DP, -vel_cap(1)/)
	  Rot_matrix(3,:)=(/ -vel_cap(2), vel_cap(1), 0.0_DP/) !!!!!!!! Do random rotationof e^_temp  around v^ by a random angle
	  Rot_angle=2.0*pi*ran2(idum)
	  Rot_matrix = Identity + sin(Rot_angle)*Rot_matrix + (1-cos(Rot_angle))*MATMUL(Rot_matrix,Rot_matrix)	
	    
	  e_cap= MATMUL(Rot_matrix,e_cap_temp) !!!!!!!!! Select the axis of rotation by rodrigues formula 
	  e_cap = e_cap/NORM2(e_cap)
	  
	   
	  
	      
	  Rot_matrix(1,:)=(/ 0.0_DP, -e_cap(3), e_cap(2)/)
	  Rot_matrix(2,:)=(/ e_cap(3), 0.0_DP, -e_cap(1)/)
	  Rot_matrix(3,:)=(/ -e_cap(2), e_cap(1), 0.0_DP/)
	  !!!! SCALAR NOISE !!!!!
	  noise_angle = eta*pi*( 2.0_DP*ran2(idum)-1.0_DP) 
	  Rot_matrix = Identity + sin(noise_angle)*Rot_matrix + (1-cos(noise_angle))*MATMUL(Rot_matrix,Rot_matrix)
	  
	  self_prop(i,:)=MATMUL(Rot_matrix,vel_cap) !!!!!!!! Do random rotation	 of v^  around e_cap by a random angle  
	  self_prop(i,:)=self_prop(i,:)/NORM2(self_prop(i,:))	
	  
	  
     !!	******************UPGRADE THE ANGLES************************
     
  

     !!************* PRINT full config. after each sample_step******
     IF(mod(step,sample_step)==0)THEN
         vx0=(position(i,1)-position_prev(i,1))/(DBLE(sample_step*h))
         vy0=(position(i,2)-position_prev(i,2))/(DBLE(sample_step*h))
         vz0=(position(i,3)-position_prev(i,3))/(DBLE(sample_step*h))
         ux0=vx0/dsqrt(vx0**2+vy0**2+vz0**2)
         uy0=vy0/dsqrt(vx0**2+vy0**2+vz0**2)
         uz0=vz0/dsqrt(vx0**2+vy0**2+vz0**2)
         
        write(21,*) sample_counter,position(i,1), position(i,2), position(i,3), ux0,uy0,uz0,i,Tg_state(i),yPZ
	  
	   !write(21,*) sample_counter,position(i,1), position(i,2), position(i,3),velocity(i,:)/NORM2(velocity(i,:)), i
	  
	  
	  
	   position_prev(i,1)= position(i,1)
	   position_prev(i,2)= position(i,2)
	   position_prev(i,3)= position(i,3)
     END IF
     !!************************************************************
     
  

        END DO !!!! END Cell Index Loop
        
  IF(mod(step,sample_step)==0)THEN   
      sample_counter=sample_counter+1
  END IF

  Ntemp=Ntot
  Volume_temp=Volume
!!!%%%%%%% INPUT OF NEW CELLS AT ADM/PNT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 IF(mod(step,step_input) == 0)THEN
	
	   do k=Ntot+1,Ntot+1 
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   
	   anglemin = (pi/2.0) - (acos(zmid/R))
	   anglemax = ( (pi/2.0) - (acos(zmid/R)) ) + 2.0*(acos(zmid/R))
	    
    60 polar_angle= anglemin +(anglemax-anglemin)*ran2(idum)   !! 30 to 150 degree angle
	   rad=(R-rc)*dsqrt(ran2(idum))
	
	   testx= rad*dcos(polar_angle)
	   testz= rad*dsin(polar_angle)
	  
	 IF (testz >= (zmid+rc) ) then 
	   position(k,1)= testx
	   position(k,3)= testz
	   position(k,2)= ymin+1.0*rc  !!!!!!  we put cells at lower layers, fixed y  
	 ELSE
	  GOTO 60
	 END IF	 
	   azimuth_angle=2.0_DP*pi*ran2(idum)
	   polar_angle= dacos(2.0*ran2(idum) - 1.0)	!!!!!!!!!! cos^-1[-1,1]   
	   self_prop(k,1)= dsin(polar_angle)*dcos(azimuth_angle)
	   self_prop(k,2)= dsin(polar_angle)*dsin(azimuth_angle)  
	   self_prop(k,3) = dcos(polar_angle)	 
	   position_prev(k,1)= position(k,1)
	   position_prev(k,2)= position(k,2)
	   position_prev(k,3)= position(k,3)
	end do
	  
   Ntot=Ntot+1
  
!!!%%%%%%%%%%%% ADJUST VOLUME TO INCORPORATE NEW CELLS %%%%%%%%%%%%
   dN=(Ntot-Ntemp)
	if(Ntot<Ntemp) dN=0
	dvol=DBLE(dN)*volume_temp/DBLE(Ntemp)
	
	dl= 2.0_DP*dvol/(pi*R**2)
	ymax=ymax+dl
	yPZ=yPZ+dl
!!!!These defines modified boundaries -- use this

volume= (0.5_DP*pi*R**2)*ymax + 0.25_DP*(4.0_DP/3.0_DP)*pi*R**3
    
END IF	
!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END DO !!!! time-loop ends



CLOSE(21)

END PROGRAM horseshoe

!!%%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















!!!!/////// CELL-CELL force subroutine////////////////////////////////////

!!! This subroutine Computes the force on a cell by other cells

subroutine cell_force(yPZ,x1,y1,z1,x2,y2,z2,r1,r2,Tg,fx0,fy0,fz0)
        USE numz
        USE constants
        IMPLICIT NONE
	
      REAL(DP),PARAMETER::del=1.2_DP
	  REAL(DP),INTENT(IN)::yPZ,x1,y1,z1,x2,y2,z2,r1,r2
	  INTEGER,INTENT(IN)::Tg
      REAL(DP),INTENT(OUT):: fx0,fy0,fz0
      REAL(DP)::x12,y12,z12,d12,req,r0,rep
      INTEGER:: ID1, ID2  !!!!!!!!!!!!!!!! ADM=1, PZ+PZ-neck=2, PSM-Right=3, PSM-Left=4
      
      fx0=0.0_DP
	  fy0=0.0_DP
	  fz0=0.0_DP
	  
	  
	  if(Tg==0)then
     rep=30.0_DP
     elseif(Tg==1)then
     rep=7.0*30_DP    !!!!! Perturbation Strength = 7 times the normal repulsion level
     end if
    
    
	  
IF( y1 <= yPZ .and. z1>= zmid)THEN
 	ID1=1
ELSEIF( y1 > yPZ )THEN
	ID1=2 
ELSEIF( y1 <= yPZ .and. z1< zmid .and. x1>0)THEN
	ID1=3
ELSEIF( y1 <= yPZ .and. z1< zmid .and. x1<0)THEN	
	ID1=4
END IF

IF( y2 <= yPZ .and. z2>= zmid)THEN
 	ID2=1
ELSEIF( y2 > yPZ )THEN
	ID2=2 
ELSEIF( y2 <= yPZ .and. z2< zmid .and. x2>0)THEN
	ID2=3
ELSEIF( y2 <= yPZ .and. z2< zmid .and. x2<r0)THEN	
	ID2=4
END IF

  
IF(ID1 == ID2 .or. ID1==2 .or. ID2==2)THEN	  
	x12=x1-x2
	y12=y1-y2
	z12=z1-z2
	d12=dsqrt(x12**2+y12**2+z12**2)
	
	req=(r1+r2)
	r0=del*req
	
	if(d12<req .and. d12 > 0.0_DP)then
	  fx0=rep*(x12/d12)*(req-d12)/req
	  fy0=rep*(y12/d12)*(req-d12)/req
	  fz0=rep*(z12/d12)*(req-d12)/req
	endif
	
	if(d12<= r0 .and. d12 >= req)then
	  fx0=-Fadh*(x12/d12)*(d12-req)/(r0-req)
	  fy0=-Fadh*(y12/d12)*(d12-req)/(r0-req)
	  fz0=-Fadh*(z12/d12)*(d12-req)/(r0-req)
	endif
	
	if(d12 > r0 .or. d12==0.0_DP)then
	  fx0=0.0_DP
	  fy0=0.0_DP
	  fz0=0.0_DP
	endif
END IF	
	return
end subroutine cell_force

!!!!/////////////////////////////////////////////////////////////////////////



!!!!/////// Reflecting boundary condition subroutine////////////////////////////////////

!!! This subroutine imposes the reflecting boundary condition on each cell

subroutine boundary(ymax,yPZ,rc,xi,yi,zi,xf,yf,zf)
        USE numz
        USE constants
        IMPLICIT NONE
	
	    REAL(DP),INTENT(IN)::ymax,yPZ,rc,xi,yi,zi
        REAL(DP),INTENT(INOUT):: xf,yf,zf
		REAL(DP):: Rf, alpha
		
         !! For ADM regions :::::::::::::::::::::::::::
IF( yi <= yPZ .and. zi>= zmid)THEN
    if(yf <= rc)  yf=2.0_DP*(rc)-yf
    if(zf <= (zmid +rc)) zf=2.0_DP*(zmid +rc)-zf
       
    Rf=dsqrt(xf**2+ zf**2)	    
	if(Rf >= (R-rc))then	    
	    alpha= atan(zf/xf)	    
	    if(alpha>0)then
	    xf= (2*(R-rc) - Rf)*dcos(alpha)
	    zf= (2*(R-rc) - Rf)*dsin(alpha) 
	    elseif(alpha<0)then
	    xf= (2*(R-rc) - Rf)*dcos(pi-abs(alpha))
	    zf= (2*(R-rc) - Rf)*dsin(pi-abs(alpha))
	    end if	   
   end if
	          
 END IF
 
  !! For PSM-Right regions :::::::::::::::::::::::::::
IF( yi <= yPZ .and. zi< zmid .and. xi>0)THEN
    if(yf <= rc)  yf=2.0_DP*(rc)-yf
    if(zf <= rc)  zf=2.0_DP*(rc)-zf   
    if(zf >= (zmid -rc)) zf=2.0_DP*(zmid -rc)-zf
    if(xf <= rc)  xf=2.0_DP*(rc)-xf
    
    Rf=dsqrt(xf**2+ zf**2)	    
	if(Rf >= (R-rc))then	    
	    alpha= atan(zf/xf)	    
	    if(alpha>0)then
	    xf= (2*(R-rc) - Rf)*dcos(alpha)
	    zf= (2*(R-rc) - Rf)*dsin(alpha) 
	    elseif(alpha<0)then
	    xf= (2*(R-rc) - Rf)*dcos(pi-abs(alpha))
	    zf= (2*(R-rc) - Rf)*dsin(pi-abs(alpha))
	    end if	   
   end if
             
 END IF
 
 !! For PSM-Left regions :::::::::::::::::::::::::::
IF( yi <= yPZ .and. zi< zmid .and. xi<0)THEN
    if(yf <= rc)  yf=2.0_DP*(rc)-yf
    if(zf <= rc)  zf=2.0_DP*(rc)-zf   
    if(zf >= (zmid-rc)) zf=2.0_DP*(zmid -rc)-zf
    if(xf >= -rc)  xf=2.0_DP*(-rc)-xf
    
    Rf=dsqrt(xf**2+ zf**2)	    
	if(Rf >= (R-rc))then	    
	    alpha= atan(zf/xf)	    
	    if(alpha>0)then
	    xf= (2*(R-rc) - Rf)*dcos(alpha)
	    zf= (2*(R-rc) - Rf)*dsin(alpha) 
	    elseif(alpha<0)then
	    xf= (2*(R-rc) - Rf)*dcos(pi-abs(alpha))
	    zf= (2*(R-rc) - Rf)*dsin(pi-abs(alpha))
	    end if	   
   end if
   
 END IF
 
 !! For ADM-PZ-neck regions :::::::::::::::::::::::::::
IF( yi > yPZ .and. yi<= ymax)THEN
    if(zf <= rc)  zf=2.0_DP*(rc)-zf
      
    Rf=dsqrt(xf**2+ zf**2)	    
	if(Rf >= (R-rc))then	    
	    alpha= atan(zf/xf)	    
	    if(alpha>0)then
	    xf= (2*(R-rc) - Rf)*dcos(alpha)
	    zf= (2*(R-rc) - Rf)*dsin(alpha) 
	    elseif(alpha<0)then
	    xf= (2*(R-rc) - Rf)*dcos(pi-abs(alpha))
	    zf= (2*(R-rc) - Rf)*dsin(pi-abs(alpha))
	    end if	   
   end if
   !! Mid-plane
   IF (zi>=(zmid -rc) .and. zi<=(zmid+rc))THEN
      if(yf<=(yPZ+rc)) yf=2.0_DP*(yPZ+rc)-yf
    END IF
   !!PSM-left-right divider
    IF (xi>=(-rc) .and. xi<=(rc) .and. zi <=(zmid +rc) )THEN
      if(yf<=(yPZ+rc)) yf=2.0_DP*(yPZ+rc)-yf
    END IF  
   
 END IF
 
 !! For DM+PZ regions :::::::::::::::::::::::::::
IF( yi > ymax )THEN
    if(zf <= rc) zf=2.0_DP*(rc)-zf
    
    Rf=dsqrt(xf**2+ (yf-ymax)**2+ zf**2)	    
	if(Rf >= (R-rc))then	    	    	    	    
	    xf= (2*(R-rc) - Rf)*(xf)/Rf
	    yf= (2*(R-rc) - Rf)*(yf -ymax)/Rf + ymax
	    zf= (2*(R-rc) - Rf)*(zf)/Rf 	       
   end if
             
 END IF

	return
end subroutine boundary

!!!!/////////////////////////////////////////////////////////////////////////




!!!!/////// Uniform Random number generators////////////////////////////////////

FUNCTION ran2(idum)
  USE numz
  IMPLICIT NONE
  REAL(DP):: ran2
  !INTEGER,INTENT(inout),OPTIONAL::idum
  INTEGER,INTENT(inout)::idum
  INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
  INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
  INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
  INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
  REAL(DP),PARAMETER::AM=1.0_DP/IM1,EPS=1.2e-7,RNMX=1.0_DP-EPS
  INTEGER::idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  IF (idum<0) THEN
     idum=MAX(-idum,1)
     idum2=idum
      DO j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         IF (idum<0) idum=idum+IM1
         IF (j.LE.NTAB) iv(j)=idum
      ENDDO
      iy=iv(1)
   ENDIF
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   IF (idum<0) idum=idum+IM1
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2
   IF (idum2<0) idum2=idum2+IM2
   j=1+iy/NDIV
   iy=iv(j)-idum2
   iv(j)=idum
   IF(iy.LT.1)iy=iy+IMM1
   ran2=MIN(AM*iy,RNMX)
   RETURN
 END FUNCTION ran2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!/////// Gaussian random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\
 Subroutine gasdev(g2,idum,mean,variance) 
  use numz
  Implicit none
      INTEGER,INTENT(INOUT)::idum
      REAL(DP),INTENT(IN)::variance,mean
      REAL(DP),INTENT(OUT)::g2
      REAL(DP)::ran2,g1
      INTEGER:: iset
      REAL(DP):: fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
    !  if (iset.eq.0) then
        DO
        v1=2.0_DP*ran2(idum)-1._DP
        v2=2._DP*ran2(idum)-1._DP
        rsq=v1**2+v2**2
        if((rsq<1.0).AND.(rsq/=0.0))EXIT
        ENDDO
        fac=variance*DSQRT(-2._DP*log(rsq)/(rsq))
        g1=v1*fac+mean
        g2=v2*fac+mean      
      END subroutine gasdev

!!!!!!!!!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
