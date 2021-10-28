C     This program computes the spanwise lift distribution and
C     induced drag of planar wings using a discrete vortex
C     representation of the wake and a concentrated bound vortex.
C     Lift and Cl distributions are written to the file Wing.dat.
C
C     © 1986  Ilan Kroo  All rights reserved.
C
C
C
C     Declarations
C     ------------

      implicit real(a-h,o-z)

C      real SWEEP, TAPER, AR, BETA, CL, CD, Cm, MACH, LSUM, MSUM, DSUM
	
	real AR, BETA, MACH, LSUM, MSUM, DSUM
      real WIDTH, X, RN, RT, RNI
      real RTI, RY, RX, RYI, WT, WTI, R0, R1, R0I, R1I, WB, WBI, Wto
      real ALPHA, W, R, RIMAGE, BARL, PI, D2R, BETAI, CLSEC, SECLFT
      integer NPANEL, i, j, j1, NOW, MAXPAN

      parameter (MAXPAN = 1000)

      REAL Y(MAXPAN), K(MAXPAN,MAXPAN), CHORD(MAXPAN), B(MAXPAN)
      REAL L(MAXPAN), X0(MAXPAN), VB(MAXPAN), etasez(MAXPAN)
      INTEGER IP(MAXPAN), READUNIT, WRITEUNIT

      PARAMETER (READUNIT = 5,  WRITEUNIT = 6)

C Luigi mod
	integer MAXSTATION,MAXALPHA,Nstation
	parameter (MAXSTATION = 100)
	parameter (MAXALPHA = 200)
	real a
	real Aper, Croot, Alfinw, Alffinw, Alfinc, Lam
        real Nalfaw,rocr
	real Xlestat(MAXSTATION), Ystat(MAXSTATION), Cord(MAXSTATION)
	real Twist(MAXSTATION), Cm02d(MAXSTATION), SWEEPcr4(MAXSTATION)
	real SWEEP(MAXSTATION), AZL(MAXSTATION), XacStat(MAXSTATION)
	real TANSW(MAXSTATION), SECEFF(MAXSTATION), COSEFF(MAXSTATION)
	real SINEFF(MAXSTATION)
	real CL(MAXALPHA), CD(MAXALPHA), Cm(MAXALPHA), E(MAXALPHA)
	real AlphaW(MAXALPHA), CLwq(MAXALPHA), ccl(MAXALPHA,MAXPAN)
	real SECL(MAXALPHA,MAXPAN), CLSTAT(MAXALPHA,MAXPAN)
	real Cmac2d(MAXPAN),A0L(MAXPAN),Xac2d(MAXPAN),Cmac2dd(MAXPAN),
     #Xac2dd(MAXPAN),xlee(MAXPAN), CHORDa(MAXPAN)
	real WSHOUT(MAXPAN),Clsez(MAXALPHA,MAXPAN)
	character*20  fileinp
	integer jj,KK
	real frecr(MAXPAN), xle(MAXPAN), Kcord, Sw, Cm0wing
	real Cladd(MAXPAN), Clbas(MAXPAN), Xac(MAXALPHA),
     #Yac(MAXALPHA),cm1(MAXALPHA),cm2(MAXALPHA),CmacW(MAXALPHA),
     #Xcp(MAXALPHA)
        real amom(MAXALPHA,MAXPAN),dta(MAXALPHA,MAXPAN)
	real damom(MAXALPHA,MAXPAN),damomdim(MAXALPHA,MAXPAN)
	real tadim(MAXALPHA,MAXPAN),amomdim(MAXALPHA,MAXPAN)
        real dtadim(MAXALPHA,MAXPAN),dMtors(MAXALPHA,MAXPAN),
     #Mtors(MAXALPHA,MAXPAN),ta(MAXALPHA,MAXPAN),tatot(MAXALPHA),
     #amomtot(MAXALPHA),qd(MAXALPHA)
	real puntap(MAXALPHA,MAXPAN),puntapd(MAXALPHA,MAXPAN)
C Luigi mod

C     Input parameters
C     ----------------
C
      PI=3.14159
      D2R=PI/180.

      WRITE(WRITEUNIT,*)'               Planar Wing Analysis Program' 
	WRITE(WRITEUNIT,*)'                     Ilan Kroo       1986'
      WRITE(WRITEUNIT,*)'              (modificato per ali generiche)'
1     WRITE(WRITEUNIT,101)
101   FORMAT(' ',/,/,/,/)
C Luigimod      lettura geometria e aerodinamica dell'ala
!	 WRITE(*,*)'Enter name file input..'
!	 READ (*,'(a)')fileinp
	 fileinp = 'weiss.inp'
	 OPEN (10,file=fileinp)
	 READ (10,'(a)') a
	 READ (10,*) Mach,AR,Croot,Wto,rocr,xtors
	 READ (10,'(a)') a
	 READ (10,*) Alfinw,Alffinw,Nalfaw
	 READ (10,'(a)') a
	 READ (10,*) NPANEL,Nstation
	 READ (10,'(a)') a
	 READ (10,*) (Ystat(i),i=1,Nstation)
	 READ (10,'(a)') a
	 READ (10,*) (Cord(i),i=1,Nstation)
	 READ (10,'(a)') a
	 READ (10,*) (Twist(i),i=1,Nstation)
	 READ (10,'(a)') a
	 READ (10,*) (AZL(i),i=1,Nstation)
	 READ (10,'(a)') a
	 READ (10,*) (XacStat(i),i=1,Nstation)
	 READ (10,'(a)') a
	 READ (10,*) (Cm02d(i),i=1,Nstation)
	 READ (10,'(a)') a
	 READ (10,*) (Xlestat(i),i=1,Nstation)
	 CLOSE(10)

!      WRITE(WRITEUNIT,*) ' Mach Number (< 1):'

!      READ(READUNIT,*) MACH
      IF(Nalfaw.EQ.1) THEN 
	    Alfinc=0
	ELSE
	    Alfinc=(Alffinw-Alfinw)/(Nalfaw-1)
	ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Stima di alcune grandezze geometriche fondamentali
	Kcord=0.		!! Superficie adimensionale !!
      DO I=1,Nstation-1
	 Kcord=Kcord+(Cord(I+1)/Croot+Cord(I)/Croot)
     #*(Ystat(I+1)-Ystat(I))
      ENDDO
	aper = AR*Kcord/2*Croot
!     Aperadi=AR*Kcord/2
	Sw=Kcord*Croot*Aper/2
!	AR = aper**2/Sw
	
	DO I=1,Nstation-1
	    DY=(Ystat(I+1)-Ystat(I))*0.5*Aper
	    DXcr4=(Xlestat(I+1)+0.25*Cord(I+1))-(Xlestat(I)+0.25*Cord(I))
	    DX=(Xlestat(I+1)-Xlestat(I))
	    SWEEPcr4(I)=ATAN(DXcr4/DY)    ! Freccia ad un quarto in [rad]
		SWEEP(I)=ATAN(DX/DY)
	ENDDO

 
      WRITE(WRITEUNIT,'(/,A,/)')' Computing...'
      
C     Calculations
C     ------------
      BETA=SQRT(1.-MACH**2)
      BETAI=1./BETA
      
	ALPHA=Alfinw
 
      WIDTH=1./NPANEL

	DO JJ=1,Nalfaw  ! Main loop over alpha!


	DO I = 1,NPANEL

          VB(I)=REAL(I)/NPANEL    !
          IF (I.EQ.1) Y(I)=VB(I)*.5               !
          IF (I.NE.1) Y(I)=(VB(I)+VB(I-1))*.5     !
		
!		DO IS=1,Nstation
!			TANSW(I) =TAN(SWEEPcr4(IS))
!		ENDDO


	    DO IS=1,Nstation-1

	        IF((Y(I).GE.Ystat(IS)).AND.(Y(I).LE.Ystat(IS+1))) THEN

	            PERC     =(Y(I)-Ystat(IS))/(Ystat(IS+1)-Ystat(IS))
	            CHORD(I) =(Cord(IS)+(Cord(IS+1)-Cord(IS))*PERC)
     #                 /(Aper/2)    !!Adimensionale rispetto semi-apertura!!
	            A0L(I)   =AZL(IS)+(AZL(IS+1)-AZL(IS))*PERC
	            WSHOUT(I)=Twist(IS)+(Twist(IS+1)-Twist(IS))*PERC
			Cmac2d(I)=Cm02d(IS)+(Cm02d(IS+1)-Cm02d(IS))*PERC
                  TANSW(I) =TAN(SWEEPcr4(IS))
		XAC2d(I) =XacStat(IS)+(XacStat(IS+1)-XacStat(IS))*PERC
			XLE(I)=(Y(I)*(Aper/2)-Ystat(IS)*(Aper/2))
     #				*TAN(SWEEP(IS))+Xlestat(IS)   !! Dimensionale!!
                  SECEFF(I)=SQRT((TANSW(I)/BETA)**2+1.)
                  COSEFF(I)=1./SECEFF(I)
                  SINEFF(I)=TANSW(I)/BETA*COSEFF(I)
	            L(I)=VB(I)*SECEFF(I)

	        END IF

          END DO

	END DO
C
C-----------------------------------------------------------------------
C	Computing influence coefficients
C
	DO 30 I=1,NPANEL
	    X=(Y(I)*TANSW(I)+CHORD(I)*.5)*BETAI
          RN=X*COSEFF(I)-Y(I)*SINEFF(I)
          RT=X*SINEFF(I)+Y(I)*COSEFF(I)
          RNI=RN+2.*Y(I)*SINEFF(I)
          RTI=RT-2.*Y(I)*COSEFF(I)
          R0=SQRT(RN*RN+RT*RT)
          R0I=SQRT(RNI*RNI+RTI*RTI)

C      	Loop over vortices
C         WB = downwash due to bound vortex
C         WT = downwash due to trailing vorticity
C         I after the variable designates image (left side) influence

	  DO 20 J=1,NPANEL
            RY=VB(J)-Y(I)
            RX=X-VB(J)*TANSW(J)*BETAI
            RYI=RY+2.*Y(I)
            WT=(1.+RX/SQRT(RX*RX+RY*RY))/RY
            WTI=(1.+RX/SQRT(RX*RX+RYI*RYI))/RYI
            R1=SQRT(RN*RN+(L(J)-RT)**2)
            WB=((L(J)-RT)/R1+RT/R0)/RN
            R1I=SQRT(RNI*RNI+(L(J)-RTI)**2)
            WBI=((L(J)-RTI)/R1I+RTI/R0I)/RNI
            K(I,J)=(WT+WB+WTI+WBI)/(4.*PI)
20	  CONTINUE
30	CONTINUE
C
C
C-----------------------------------------------------------------------
C       L-U Decomposition of the linear system
	CALL DECOMP(NPANEL,MAXPAN,K,IP)
C-----------------------------------------------------------------------
C
C Luigi mod
	 DO I=1,NPANEL
		B(I) = (ALPHA -A0L(I) + WSHOUT(I))*D2R
	 ENDDO
C Luigi mod
C
C-----------------------------------------------------------------------
C       Back substitution for solution to system 
	CALL SOLVE(NPANEL,MAXPAN,K,B,IP)
C-----------------------------------------------------------------------
C
C	Computing drag and section loading
C
	DO 50 I=1,NPANEL
           LSUM=0
           DO 48 J=1,NPANEL+1-I
              J1=NPANEL+1-J
              LSUM=LSUM+B(J1)
48         CONTINUE
           SECL(JJ,I)=2.*LSUM
50      CONTINUE

C       Integrations for CL,Cm,CD
        LSUM=0
        MSUM=0
        DSUM=0
        DO 60 I=1,NPANEL
           W=0
           DO 55 J=1,NPANEL
C            The downwash produced by trailing vorticity is:
             R=(VB(J)-Y(I))*4.*PI
             RIMAGE=(VB(J)+Y(I))*4.*PI
             W=W+B(J)*(1./R+1./RIMAGE)
55         CONTINUE
           BARL=SECL(JJ,I)
           LSUM=LSUM+BARL*WIDTH
           MSUM=MSUM-BARL*WIDTH*Y(I)*TANSW(I)
           DSUM=DSUM+BARL*WIDTH*W
	     CLSTAT(JJ,I)=SECL(JJ,I)/CHORD(I)
60      CONTINUE
        
        CD(JJ)=DSUM*AR*.5
        CL(JJ)=LSUM*AR*.5
        CM(JJ)=MSUM*AR*AR*.25
        E(JJ) = 0.
        IF(CD(JJ).NE.0.)E(JJ)=CL(JJ)*CL(JJ)/(PI*AR*CD(JJ))

      
	AlphaW(JJ)=ALPHA
	ALPHA=ALPHA+Alfinc
	ENDDO		!!!Fine main loop over alpha!!!

      Lam=(Cord(Nstation)*Aper/2)/(Cord(1)*Aper/2)	!! Rastremazione !!

C Stima carico basico e carico addizionale
C	Queste funzioni le si utilizzano con il programma modificato 
C	di Kroo per valutare il carico basico, il carico addizionale
C	e tutte le altre grandezze che possono essere valutate sempre
C	in campo lineare (moto a potenziale ==> solo pre-stallo)
!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Nel codice di Kroo le lunghezze sono tutte adimensionalizzate 
C rispetto la semi-apertura alare; di seguito al contrario si 
C preferisce adimensionalizzare Corde e lunghezze lungo 'x' con la
C corda di radice.
	DO I=1,NPANEL
		CHORD(I)=CHORD(I)*(Aper/2)/Croot
		XLE(I)=XLE(I)/Croot
	ENDDO

cfab09-1
	CHORD(0)=cord(1)/Croot
	CHORD(NPANEL+1)=cord(NStation)/Croot

	XLE(0)=xlestat(1)/Croot
	XLE(NPANEL+1)=xlestat(Nstation)/Croot

	Y(0)=0
	Y(NPANEL+1)=1.0	

	Xac2d(0)=Xacstat(1)
	Xac2d(NPANEL+1)=Xacstat(Nstation)

	Cmac2d(0)=Cm02d(1)
	Cmac2d(NPANEL+1)=Cm02d(Nstation)
cfab09-1


C Stima corda media
	sumc=0
	sumxc=0
	sumyc=0
	
cfab09-2
	sumc=sumc+0.5*(CHORD(0)**2+CHORD(1)**2)*Y(1)
	sumxc=sumxc+0.5*(CHORD(0)*xle(0)+CHORD(1)*xle(1))*Y(1)
	sumyc=sumyc+0.5*(CHORD(0)*Y(0)+CHORD(1)*Y(1))*Y(1)
cfab09-2

	DO i=1,NPANEL-1
		dy=(Y(i+1)-Y(i))
		sumc=sumc+.5*(CHORD(i+1)**2+CHORD(i)**2)*dy  
C	Per stimare X_cmed
		sumxc=sumxc+.5*(CHORD(i+1)*xle(i+1)+CHORD(i)*xle(i))*
     #	dy
C	Per stmare Y_cmed	
		sumyc=sumyc+.5*(CHORD(i+1)*Y(i+1)+CHORD(i)*Y(i))*dy
	ENDDO

cfab09
	dy=Y(NPANEL+1)-Y(NPANEL)
	sumc=sumc+0.5*(CHORD(NPANEL)**2+CHORD(NPANEL+1)**2)*dy
	delta_xc=(CHORD(NPANEL)*XLE(NPANEL)+CHORD(NPANEL+1)*XLE(NPANEL+1))
	delta_xc=delta_xc*0.5*dy
	sumxc=sumxc+delta_xc
	delta_yc=(CHORD(NPANEL)*Y(NPANEL)+CHORD(NPANEL+1)*Y(NPANEL+1))
	sumyc=sumyc+0.5*delta_yc*dy
cfab09

	Cmed=2/kcord*sumc
	Xcmed=2/kcord*sumxc
	Ycmed=2/kcord*sumyc
	WRITE(*,*)"Cmed/Cr=",Cmed
	WRITE(*,*)"X-Cmed/Cr= ",Xcmed
	WRITE(*,*)"Y-Cmed/(b/2)= ",Ycmed
!	PAUSE
	IF (Nalfaw.GE.2.) THEN
	DO i=1,NPANEL
		Cladd(I)=(CLSTAT(2,I)-CLSTAT(1,I))/(CL(2)-CL(1))
		Clbas(I)=CLSTAT(1,I)-Cladd(I)*CL(1)
	ENDDO				
C Per avere il carico addizionale effettivo bisogna fare il prodotto
C Cl_a=Cl_add*CL (vedi anche metodo Anderson)  		
	Cm0wing=0.

cfab09
c ESTRAPOLAZIONE radice e tip
	Cladd(0)=Cladd(1)+(Cladd(1)-Cladd(2))/(Y(2)-Y(1))*Y(1)
	Clbas(0)=Clbas(1)+(Clbas(1)-Clbas(2))/(Y(2)-Y(1))*Y(1)

	dytip=Y(NPANEL)-Y(NPANEL-1)
	dyest=1.0-Y(NPANEL)
	dcladd=Cladd(NPANEL)-Cladd(NPANEL-1)
	dclbas=Clbas(NPANEL)-Clbas(NPANEL-1)
	Cladd(NPANEL+1)=Cladd(NPANEL)+dCladd/dytip*dyest
	Clbas(NPANEL+1)=Clbas(NPANEL)+dClbas/dytip*dyest
cfab09

cfab09
	DO jj=1,Nalfaw
	dyr=(Y(2)-Y(1))
	CLSTAT(jj,0)=CLSTAT(jj,1)+(CLSTAT(jj,1)-CLSTAT(jj,2))/dyr*Y(1)
	SECL(jj,0)=SECL(jj,1)+(SECL(jj,1)-SECL(jj,2))/dyr*Y(1)

	dytip=Y(NPANEL)-Y(NPANEL-1)
	dyest=1.0-Y(NPANEL)
	dclstat=CLSTAT(jj,NPANEL)-CLSTAT(jj,NPANEL-1)
	dsecl=SECL(jj,NPANEL)-SECL(jj,NPANEL-1)

	CLSTAT(jj,NPANEL+1)=CLSTAT(jj,NPANEL)+dclstat/dytip*dyest
	SECL(jj,NPANEL+1)=SECL(jj,NPANEL)+dsecl/dytip*dyest

	ENDDO				
cfab09



     	DO jj=1,Nalfaw
		sum=0
		sumcm1=0
		sumcm2=0
		sumclad=0
		sumcla=0
		sumXac=0
		sumYac=0

cfab09		DO I=1,NPANEL-1 
		DO I=0,NPANEL
		cladi1=Cladd(I+1)*CL(JJ)
		cladi=Cladd(I)*CL(JJ)
		dy=(Y(I+1)-Y(I))
C     Per stimare CmacTot e Xcp
	sumcm1=sumcm1+.5*(-Clbas(I+1)*CHORD(I+1)*(Xac2d(I+1)*CHORD(I+1)
     #	+xle(I+1))-Clbas(I)*CHORD(I)*(Xac2d(I)*CHORD(I)+xle(I)))*dy
		sumcm2=sumcm2+.5*(cmac2d(I+1)*CHORD(I+1)**2+cmac2d(I)
     #	*CHORD(I)**2)*dy
C     Per stimare Xac
		sumXac=sumXac+.5*(cladi1*CHORD(I+1)*(Xac2d(I+1)*CHORD(I+1)
     #	+xle(I+1))+cladi*CHORD(I)*(Xac2d(I)*CHORD(I)+xle(I)))*dy
C     Per stimare Yac
	sumYac=sumYac+.5*(cladi1*CHORD(I+1)*Y(i+1)+cladi*CHORD(I)
     #	*Y(I))*dy
   		ENDDO
CCCC  Fine calcolo integrali e valutazione grandezze globali
		IF(CL(jj).eq.0.) GOTO 1000
  		Xac(jj)=2/(CL(jj)*kcord)*sumXac*Croot
		Yac(jj)=2/(CL(jj)*kcord)*sumYac*Aper/2
		cm1(jj)=2/(kcord*cmed)*sumcm1
		cm2(jj)=2/(kcord*cmed)*sumcm2
		CmacW(jj)=cm1(jj)+cm2(jj)
		Xac_cr=Xac(jj)
		Yac_bm=Yac(jj)
	    Xac_mac=(Xac_cr-Xcmed*Croot)/(Cmed*Croot)
1000	ENDDO

	ELSE
	WRITE(*,*) "Impossibile valutare carico basico e addizionale..."
      WRITE(*,*) "Inserire almeno due valori di Alpha..."
!	PAUSE
	ENDIF

C Sezione file di output

	OPEN(14,file='WingSTAT.out')
	WRITE (14,37)'%Alpha[deg]','c[m]','y/(b/2)','Xle/Cr','Clstaz',
     #'Cl_bas','Cl_add','Ccl/cavg','Cmac2d'


cfab09	DO I=1,NPANEL
	DO I=0,NPANEL+1
		  CHORD(I)=CHORD(I)*Croot	!! Corda dimensionale !!
	ENDDO
	DO JJ=1,Nalfaw
	   DO I=0,NPANEL+1
	      WRITE (14,'(9(f10.5))') AlphaW(JJ),CHORD(I),Y(I),xle(I)
     #,CLSTAT(JJ,I),Clbas(I),Cladd(I)*CL(JJ),SECL(JJ,I)*AR*0.5,cmac2d(I)
	   ENDDO
	WRITE (14,'(a)')
      WRITE (14,'(a)')
	ENDDO
	CLOSE(14)

	OPEN (unit=15,file='weiss.out')

35     FORMAT(500(f12.6))
36	FORMAT(4x,a,9x,a,9x,a,7x,a,7x,a,10x,a)
37	format(a,3x,a,5x,a,3x,a,4x,a,4x,a,4x,a,4x,a,4x,a)
911   format(a1,a28,F12.6)

      write(15,911) '%','                 AR=',AR
      write(15,911) '%','                  b=',Aper
      write(15,911) '%','                 Sw=',Sw
      write(15,911) '%','                 C_radice=',Croot
      write(15,911) '%','Rapp rastremazione =',Lam

      write(15,'(A1,A20)') '%','   ****   MAC   ****    '     
      write(15,911) '%','                MAC/Cr=',cmed
      write(15,911) '%','                MAC [m]=',cmed*Croot
      write(15,911) '%','              X_MAC [m]=',xcmed*Croot
      write(15,911) '%','              Y_MAC [m]=',ycmed*Aper/2

      write(15,'(A1,A25)') '%','  ****  AER. CENTER   ****  '     
      write(15,911) '%','POS in x         X_AC [m]=',Xac_cr
      write(15,911) '%','POS in aper. y   Y_AC [m]=',Yac_bm
	write(15,911) '%','POS in x (%MAC)  X_AC/MAC=',Xac_mac

      write(15,'(A1,A25)') '%','  ****  WING AC MOM COEFF   ****  '     
      write(15,911) '%','basic load contr  Cmac1=',cm1(JJ-1)
      write(15,911) '%','local cmac contr  Cmac2=',cm2(JJ-1)
      write(15,911) '%','        WING      CMACw=',CmacW(JJ-1)

	WRITE(15,36) '%Alpha','CL','CDi','Cm(c/4)','e'
      DO JJ=1,Nalfaw
		write(15,35) AlphaW(JJ),CL(JJ),CD(JJ),CM(JJ),E(JJ)
      ENDDO
	CLOSE(15)
C!!!!!!!!!!!!!!!!!!!WEISS MODIFICA PADUANO!!!!!!!!!!!!!!!!!!
C      Calcolo Taglio e Momento
        Do i=1,nalfaw
        do j=0,npanel+1
        ccl(i,j)=SECL(i,npanel-j+1)*AR*0.5*cmed*Croot
        etasez(j)=Y(npanel-j+1)
        CLwq(i)=CL(i)
        Clsez(i,j)=CLSTAT(i,npanel-j+1)
        CHORDa(j)=CHORD(npanel-j+1)
        xlee(j)=xle(npanel-j+1)*Croot
        xac2dd(j)=xac2d(npanel-j+1)*Chorda(j)
        cmac2dd(j)=cmac2d(npanel-j+1)

	if(CL(i).eq.0) then
	CLwq(i)=0.1
	end if

	if(CL(i).lt.0) then
	CLwq(i)=-CL(i)
	end if
        qd(i)=Wto/(Sw*CLwq(i))
	Vcr=(2*qd(i)*9.81/rocr)**0.5
        end do
        end do

	ta(1,1)=0.0
	tadim(1,1)=0.
	tatot(1)=0.0
	dtadim(1,1)=0.
	dta(1,1)=0.

	DO j=1,Nalfaw
           DO i=1,NPANEL+1

                dtadim(j,i)=0.25*Aper*qd(j)*(etasez(i-1)-etasez(i))
     #          *(ccl(j,i-1)+ccl(j,i))
	        tadim(j,i)=tadim(j,i-1)+dtadim(j,i)

	        dta(j,i)=dtadim(j,i)/(qd(j)*Sw*0.5)
	        ta(j,i)=ta(j,i-1)+dta(j,i)

                tatot(j)=tatot(j)+dta(j,i)
           END DO
	END DO

	amom(1,1)=0.0
	amomtot(1)=0.0
	amomdim(1,1)=0.0
	
	DO j=1,Nalfaw
	DO i=1,NPANEL+1

      damomdim(j,i)=0.25*Aper*(etasez(i-1)-etasez(i))
     #*(tadim(j,i-1)+tadim(j,i))
	amomdim(j,i)=amomdim(j,i-1)+damomdim(j,i)

	damom(j,i)=damomdim(j,i)/(qd(j)*Sw*Aper*0.25)
	amom(j,i)=amom(j,i-1)+damom(j,i)
        amomtot(j)=amomtot(j)+damom(j,i)


        END DO
	END DO


	OPEN(19,File='wingint.res')

        write(19,*) '# GRANDEZZE ADIMENSIONALI'
	write(19,*)
	write(19,*) '# Cl        -> adim. risp. a Q e Sw'
	write(19,*) '# dT e T    -> adim. risp. a Q e Sw/2'
	write(19,*) '# dM e M    -> adim. risp. a Q , Sw/2 e b/2'
	write(19,*)
	write(19,*) '# Taglio totale adim  =',tatot(1)         !!mi interessa solo il primo angolo!!
	write(19,*) '# Momento totale adim =',amomtot(1)       !!ADAS usa da questo programma solo il primo!!
        write(19,*)

c	write(19,110)
	write(19,'(1x,A89)') '#    eta          c[m]          Cl           cCl
     >  dT           T            dM            M'

        j=1
	DO i=0,NPANEL+1
	write(19,910) etasez(i),CHORDa(i),Clsez(j,i),ccl(j,i),
     #dta(j,i),ta(j,i),damom(j,i),amom(j,i)
     
	END DO

	CLOSE(19)

	OPEN(17,File='wing-tm.out')

	write(17,*) '# RISULTATI  DIMENSIONALI  '
	write(17,*)
      write(17,*) '# Peso                  [Kg]    =', Wto
      write(17,*) '# ro quota considerata  [Kg/m3] =', rocr
      write(17,*) '# Velocita              [m/s]   =', Vcr
	write(17,*) '# pressione dinamica    [kg/mq] =', qd(1)   !!mi interessa solo il primo angolo!!
	write(17,*) '# semi apertura (b/2)   [m]     =', 0.5*Aper
	write(17,*)

      write(17,'(1x,A101)') '#    y[m]       c[m]       Cl         c*Cl      T[kg]
     >    M[kg*m]   PA_perc   PA_dim[m]   dMt[kg*m]  Mtors[kg*m]'


	MTORS(1,1)=0
	dMtors(1,1)=0
	puntap(1,1)=999
	puntapd(1,1)=999
	
	DO j=1,Nalfaw
	DO i=1,NPANEL+1


	if(Clsez(j,i).gt.0) then
	abcl=Clsez(j,i)
	else
	abcl=-Clsez(j,i)
	end if


	if(abcl.gt.0.005) then
	puntap(j,i)=(Xac2dd(i))/CHORDa(i)  !!xac definito in weiss Š gi… privato di Xle rispetto a multhopp!!
     #  *100.-cmac2dd(i)/Clsez(j,i)*100.
	puntapd(j,i)=xlee(i)+CHORDa(i)*puntap(j,i)/100.
	else
	puntap(j,i)=999
	puntapd(j,i)=999
	end if

	if(i.eq.2) then
	pa=puntapd(j,i)
	else
	pa=(puntapd(j,i)+puntapd(j,i-1))*0.5
	end if
	dMtORS(j,i)=dtadim(j,i)*(xtors-pa)
	Mtors(j,i)=Mtors(j,i-1)+dMtors(j,i)
	END DO
        END DO
        
        j=1
	DO i=0,NPANEL+1
        ydim=etasez(i)*Aper*0.5
	ccls=Clsez(j,i)*CHORDa(i)
	write(17,931) ydim,CHORDa(i),Clsez(j,i),ccls,tadim(j,i),
     > amomdim(j,i),puntap(j,i),puntapd(j,i),dMtors(j,i),Mtors(j,i)
!!     >,xac2dd(i),xlee(i)
	END DO


	CLOSE(17)

c 100   format(113H #    y[m]       c[m]       Cl         c*Cl      T[kg]
c     +    M[kg*m]   PA_perc   PA_dim[m]   dMt[kg*m]  Mtors[kg*m])

c 110   format(99H #    eta          c[m]          Cl           cCl
c     +  dT           T            dM            M)

910   format(1x,8(F10.5,3x))
931   format(1x,14(F14.3,2x))
	END


C-----------------------------------------------------------------------
C	MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
	SUBROUTINE DECOMP(N,NDIM,C,IP)
	REAL C(NDIM,NDIM),T
	INTEGER IP(NDIM)
	IP(N)=1
	DO 6 K=1,N
	IF(K.EQ.N) GOTO 5
	KP1=K+1
	M=K
	DO 1 I=KP1,N
	IF(ABS(C(I,K)).GT. ABS(C(M,K)))M=I
1	CONTINUE
	IP(K)=M
	IF(M.NE.K)IP(N)=-IP(N)
	T=C(M,K)
	C(M,K)=C(K,K)
	C(K,K)=T
	IF (T.EQ.0.0)GOTO 5
	DO 2 I=KP1,N
	C(I,K)=-C(I,K)/T
2	CONTINUE
	DO 4 J=KP1,N
	T=C(M,J)
	C(M,J)=C(K,J)
	C(K,J)=T
	IF (T.EQ.0.0)GOTO 4
	DO 3 I=KP1,N
	C(I,J)=C(I,J)+C(I,K)*T
3	CONTINUE
4	CONTINUE
5	IF(C(K,K).EQ.0.0)IP(N)=0
6	CONTINUE
	RETURN
	END

C-----------------------------------------------------------------------
C	SOLUTION OF LINEAR SYSTEM CX = B
	SUBROUTINE SOLVE(N,NDIM,C,B,IP)
	REAL C(NDIM,NDIM),B(NDIM),T
	INTEGER IP(NDIM)
	IF(N.EQ.1)GOTO 9
	NM1=N-1
	DO 7 K = 1,NM1
	KP1=K+1
	M=IP(K)
	T=B(M)
	B(M)=B(K)
	B(K)=T
	DO 7 I=KP1,N
	B(I)=B(I)+C(I,K)*T
7	CONTINUE
	DO 8 KB=1,NM1
	KM1=N-KB
	K=KM1+1
	B(K)=B(K)/C(K,K)
	T=-B(K)
	DO 8 I=1,KM1
	B(I)=B(I)+C(I,K)*T
8	CONTINUE
9	B(1)=B(1)/C(1,1)
10	RETURN
	END

