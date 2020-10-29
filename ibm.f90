! program ibm.f90

! reproduction gene by gene

! Uses periodic boundary conditions
! Number of individuals grows until N(t) >= N

! let S increases up to S + nrmax to find surrugate mother
! if no new mother is available within S+nrmax, individual dies.

! identical to topobar6 but also
! saves information on time to common ancestor for phylogeny and genealogy

! Rectangular area

! Marcus A.M. de Aguiar  - 01/Jul/2016 -

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Reproduction 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Each individual can have one or two offspring
!
! Conditions for individual k have 2 offspring:
! 1) total number of individuals < nco
! 2) number of indiv. in pi*S^2 around k < (average-density)*(free sites around indiv.)
!    Obs. if the number of free sites around k < pi*S^2 (too many forbiden sites) area 
!         is increased to get a better estimate of the local density
! 3) mnpm >= P
! 4) number of individuals at position x(k),y(k) < 5
!
! When these conditions are satisfied k chooses a kmate among its compatible mates (there
! are at least P of them) and has two offspring, varying the position of the recombination.
!
! Conditions for individual k have 1 offspring
! If the conditions for 2 offspring fail, then we check if:
! 1) aux > Q   - probability Q of not leaving descendants
! 2) mnpm >= P - at least P compatible mates in S
! In this case a kmate is chosen among the compatible neighbors of k.
! Otherwise another 'mother' in the neighborhood of k is chosen to be reproduce in its place.
!
! If another 'mother' cannot be found, k dies without being replaced and the population decreases.
! Conditions of new mother:
! 1) k must have at least two spatial neighbors in S or in S+nrmax
! 2) new mother must have at least one compatible neighbor in its vicinity S+nrmax
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! module defining global variables
MODULE globals
INTEGER(1), ALLOCATABLE, SAVE :: g(:,:)
INTEGER(2), SAVE :: iibound
INTEGER(4), SAVE :: ndim,nf1,nf2,nb,mnpm,ineighborg,ineighbor,nrmax,igt,iitime
INTEGER(4), SAVE :: nc,nct,nco,ineighborden
INTEGER(2), ALLOCATABLE, SAVE :: x(:),y(:),neig(:),neigsp(:),ispecies(:,:)
INTEGER(2), ALLOCATABLE, SAVE :: worg(:,:,:),nworg(:,:)
INTEGER(4), ALLOCATABLE, SAVE :: ndq(:),idqx(:,:),idqy(:,:),ispv(:)
REAL, SAVE :: rg,qmat,aux,rho0,radius,mut,diff
INTEGER(4), ALLOCATABLE, SAVE :: t(:,:),tp(:,:),p1(:),p2(:),ispidx(:)  ! linha nova para filogenia
INTEGER(4), ALLOCATABLE, SAVE :: nwr   ! numero de individuos para genealogia
END MODULE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module that reads and writes the populations
MODULE readwrite
INTEGER, SAVE :: iread
CHARACTER*25, SAVE :: name
CHARACTER*30, SAVE :: name1
CONTAINS
    
	SUBROUTINE READPOP
	USE globals
	! initialize or read population
	IF(iread==1) THEN
		!WRITE(6,*) ' '
		!WRITE(6,*) 'ENTER CODE NAME FOR POP FILE'    !continue a previous run, from t=iitime
		!READ(5,*) name
        	name1 = 'pop.dat'
		OPEN(UNIT=10,FILE=name1,STATUS='old')
		READ(10,*) iitime,nct,nf1,nf2
		READ(10,*) amutd,diffd !dumb only
		READ(10,*) radiusd,rgd,nbd,mnpmd !dumb only 
		DO i=1,nct
			READ (10,901) x(i),y(i),(g(i,j),j=1,nb)
		END DO
       		CLOSE(10)
        	name1 = 'mrcat.dat'  !!write file with time matrix
        	OPEN(UNIT=10,FILE=name1,STATUS='old')
        	READ(10,*) nct
        	DO i=1,nct-1
            		READ(10,904) (t(i,jj),jj=i+1,nct)
        	END DO
       		t = t + transpose(t)
        	CLOSE(10)
    	ELSE
        	nct = nco
        	DO i=1,nct   ! random distribution
            		ictr = 0
            		DO WHILE (ictr == 0)
                		CALL RANDOM_NUMBER(aux)
                		ii = int(aux*nf1)+1
                		CALL RANDOM_NUMBER(aux)
                		jj = int(aux*nf2)+1
                		x(i) = ii
                		y(i) = jj
                		ictr = 1
            		END DO
        	END DO
	END IF
	iitime = 0
	901 FORMAT(i4,1x,i4,1x,200000i1)
	904 FORMAT(5000(i5,1x)) !time matrix
	END SUBROUTINE READPOP

	SUBROUTINE WRITEPOP
	USE globals
	!WRITE(6,*) ' '  
	!WRITE(6,*) 'ENTER CODE NAME FOR OUTPUT POP FILE'
	!WRITE(6,*) '         '
	!READ(5,*) name
	name1 = 'pop.dat'

	OPEN(UNIT=10,FILE=name1,STATUS='UNKNOWN',POSITION='REWIND')
	WRITE(10,*) iitime,nct,nf1,nf2
	WRITE(10,*) mut,diff
	WRITE(10,*) radius,rg,nb,mnpm
	DO i=1,nct
		WRITE (10,901) x(i),y(i),(g(i,j),j=1,nb)
	END DO
	CLOSE(10)

    	name1 = 'mrcat.dat'      !!write file with time matrix
    	OPEN(UNIT=9,FILE=name1,STATUS='unknown')
    	WRITE(9,*) nct
    	DO i=1,nct-1
        	WRITE(9,904) (t(i,jj),jj=i+1,nct)
    	END DO
    	CLOSE(9)

	901 FORMAT(i4,1x,i4,1x,200000i1)
	904 FORMAT(5000(i5,1x))
	END SUBROUTINE WRITEPOP


    	SUBROUTINE WRITEPHYLOGENY
    	USE globals

    	OPEN(unit=15,file='mrcat-phy.dat',status='unknown')
    	WRITE(15,*) igt,iitime,nb
    	DO ii=1,igt
        	WRITE(15,902) (t(ispecies(ii,1),ispecies(jj,1)),jj=1,igt)
    	END DO
    	CLOSE(15)


    	OPEN(unit=15,file='phy-genetics.dat',status='unknown')
    	DO ii=1,igt
        	WRITE (15,903) (g(ispecies(ii,1),jj),jj=1,nb)
    	END DO
    	CLOSE(15)

	902 FORMAT(1000(1x,i5))
	903 FORMAT(200000i1)
    	END SUBROUTINE WRITEPHYLOGENY

END MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! begin main program
PROGRAM ibm
USE globals
USE readwrite

INTEGER(1), ALLOCATABLE :: gl(:,:)
INTEGER(2), ALLOCATABLE :: indxorg(:),xl(:),yl(:)
INTEGER(4) deltat
INTEGER(2) xj(32),yj(32)
INTEGER(2) iix,iiy,ix1,iy1,idensity

! Read input data
! Note: nf1 >= nf2
OPEN(UNIT=7,FILE='input.in',STATUS='OLD',POSITION='REWIND')
READ(7,*) ntime,nco,nf1,nf2
READ(7,*) mut,diff,qmat,deltat
READ(7,*) radius,rg,nb,mnpm
READ(7,*) iread,njump
READ(7,*) nrmax
CLOSE(7)

nrmax = nrmax + 1

! diffusion table
CALL JUMPTABLE(xj,yj)

! initialize random number generator
CALL init_random_seed()

! open species plot abundance output files
OPEN(UNIT=8,FILE='speciesplot.dat',STATUS='unknown')
OPEN(UNIT=15,FILE='species-time.dat',STATUS='unknown')


! model parameters and initializations
idensity = 100  ! maximum number of individuals per site before
     		! mating is idensity/3; after mating, and before diffusion, it might be larger
rho = float(nco)/float(nf1*nf2)
rho0 = 0.83*rho ! expected density of individuals
nc = 1.05*nco  !allocate vectors with a margin for nct > nco

!
! the factor 0.83 = 5/6 compensates for the 0.6 in the local
! density for the second offspring -- that should be 0.5


ALLOCATE (g(nc,nb),gl(nc,nb))
ALLOCATE (x(nc),y(nc),xl(nc),yl(nc))
ALLOCATE (neig(nc),neigsp(nc))
ALLOCATE (ispv(nc),ispecies(nc,nc))
ALLOCATE (worg(nf1,nf2,idensity),nworg(nf1,nf2),indxorg(nc))
ALLOCATE (t(nc,nc),tp(nc,nc),p1(nc),p2(nc),ispidx(nc))       ! linha nova para filogenia

g = 0  ! identical genomes
gl = 0
x = 0
y = 0
xl = 0
yl = 0
t = 0  ! initial times to common ancestor
tp = 0 ! aux variables
p1 = 0 ! first parent
p2 = 0 ! second parent

! initialize populations
CALL READPOP

WRITE(6,*) 'total area =',nf1*nf2
WRITE(6,*) 'average density =', rho0
WRITE(6,*) 'average number of individuals in S =',3.1416*rho*radius**2
WRITE(6,*) 'initial number of individuals =',nct

! Place the organisms by location in the poster space
nworg =0
worg = 0
indxorg = 0
DO i=1,nct     
	nworg(x(i),y(i)) = nworg(x(i),y(i)) + 1
	worg(x(i),y(i),nworg(x(i),y(i))) = i
	indxorg(i) = nworg(x(i),y(i))
END DO

! Set up template for the locations that are part of the neighborhood for a radius
ndim = 4*(radius+nrmax)**2
ALLOCATE (idqx(nrmax,ndim),idqy(nrmax,ndim),ndq(nrmax))
idqx = 0
idqy = 0
ndq = 0
DO irad = 1,nrmax 
	idq = 0
	rs = radius + irad - 1
	rs2 = rs**2
	nj = int(rs)+1
	DO itx=-nj,nj               
		DO ity=-nj,nj      
			IF(itx*itx+ity*ity <= rs2) THEN
				idq = idq + 1
				idqx(irad,idq) = itx
				idqy(irad,idq) = ity
			END IF
		END DO
	END DO
	ndq(irad) = idq
END DO

write(6,*) 'initial pop=',nct
write(6,*)
write(6,*) 'partial time   total time   total pop.  species   ind./species' 
write(6,*)

!
! Time evolution: mating, mutation and diffusion
!
ntest = nco
DO j=1,ntime
	!write(6,*) 't =',j
	!Mating
    	knext = 0   ! count individuals of the next generation
	looppop: DO k=1,nct
			CALL FINDNEIG(k,1)
        		nc0k = rho0*iibound*0.6
     			! check if extra offspring is produced
			IF(nct < ntest .AND. ineighbor < nc0k) THEN
				nthere = nworg(x(k),y(k))
				IF(ineighborg >= mnpm .AND. nthere < idensity/3) THEN
					knext = knext + 1
     					CALL RANDOM_NUMBER(aux)
     					kmate = neig(int(aux*ineighborg)+1)
               				DO kc=1,nb
               					CALL RANDOM_NUMBER(aux)
               					IF(aux < 0.5) THEN
                      					gl(knext,kc) = g(k,kc)
               					ELSE
                       					gl(knext,kc) = g(kmate,kc)
               					END IF
               					CALL RANDOM_NUMBER(aux)
               					IF(aux < mut) gl(knext,kc) = 1-gl(knext,kc)   ! Mutation
               				END DO
              				xl(knext) = x(k)
					yl(knext) = y(k)
               				! save parents
               				p1(knext) = k
               				p2(knext) = kmate
					! regular offspring - same pair k-kmate 
               				knext = knext + 1
               				DO kc=1,nb
               					CALL RANDOM_NUMBER(aux)
               					IF(aux < 0.5) THEN
                       					gl(knext,kc) = g(k,kc)
               					ELSE
                       					gl(knext,kc) = g(kmate,kc)
               					END IF
               					CALL RANDOM_NUMBER(aux)
               					IF(aux < mut) gl(knext,kc) = 1-gl(knext,kc)   ! Mutation
              				END DO
               				xl(knext) = x(k)
               				yl(knext) = y(k)
               				! save parents
               				p1(knext) = k
               				p2(knext) = kmate
				END IF
			ELSE
				! only regular offspring
				irad = 1
				DO WHILE(ineighbor < 2)
					irad = irad + 1
					IF (irad == nrmax + 1) THEN
               					CYCLE looppop
					END IF
					CALL FINDNEIG(k,irad)
				END DO
				kmother = k 
				irad = 1
				CALL FINDMATE(k,kmother,kmate,irad)
				IF(ineighborg < 2) THEN
					CYCLE looppop
				ELSE
               				knext = knext + 1
               				DO kc=1,nb
              					CALL RANDOM_NUMBER(aux)
               					IF(aux < 0.5) THEN
                       					gl(knext,kc) = g(kmother,kc)
               					ELSE
                       					gl(knext,kc) = g(kmate,kc)
               					END IF
               					CALL RANDOM_NUMBER(aux)
               					IF(aux < mut) gl(knext,kc) = 1-gl(knext,kc)   ! Mutation
               				END DO
               				xl(knext) = x(k)
               				yl(knext) = y(k)
               				! save parents
               				p1(knext) = kmother
               				p2(knext) = kmate
				END IF
			END IF
	END DO looppop

	!update number of ind. and genomes
	nct = knext
	g = gl
	gl = 0
	
	! update time matrix
	DO ii=1,nct-1
		DO jj=ii+1,nct
			it11 = t(p1(ii),p1(jj))
			it12 = t(p1(ii),p2(jj))
			it21 = t(p2(ii),p1(jj))
			it22 = t(p2(ii),p2(jj))
			itlm = min(it11,it12,it21,it22)
			tp(ii,jj) = itlm + 1
			tp(jj,ii) = itlm + 1
		END DO
	END DO
	t = tp
	tp = 0

	!Diffusion
	IF(diff /= 0.0) THEN
		DO i=1,nct 
			CALL RANDOM_NUMBER(aux)
			IF(aux < diff) THEN
				ibound = 0	
				DO WHILE(ibound == 0)
					CALL RANDOM_NUMBER(aux)
					jjump = INT(njump*aux)+1
					iix = xl(i)+xj(jjump)
					iiy = yl(i)+yj(jjump)
					IF (iix > nf1) iix = iix - nf1
					IF (iix < 1) iix = nf1 + iix
					IF (iiy > nf2) iiy = iiy - nf2
					IF (iiy < 1) iiy = nf2 + iiy
					xl(i) = iix
					yl(i) = iiy
					ibound = 1
				END DO
			END IF	
		END DO
	END IF

	x = xl
	y = yl


	! Place the organisms by location in the poster space
	nworg =0
	worg = 0
	indxorg = 0
	DO i=1,nct
		nworg(x(i),y(i)) = nworg(x(i),y(i)) + 1
		worg(x(i),y(i),nworg(x(i),y(i))) = i
		indxorg(i) = nworg(x(i),y(i))
	END DO

	! calculate species every deltat starting at jtsample
	teste = (float(j)/float(deltat) - j/deltat)
	IF (teste == 0.0 .AND. j >= deltat) THEN
        	CALL FINDSPECIES
        	write(15,*) j+iitime, igt
        	write(6,*) j,j+iitime,nct,igt,(ispv(k),k=1,igt)
	END IF

! loop in time
END DO
CLOSE(15)
CLOSE(19)

DO i=1,igt
	ij = ispv(i)
	DO jjj=1,ij
		jj = ispecies(i,jjj)
		WRITE(8,*) x(jj),y(jj),i
	END DO
END DO

OPEN(unit=19,file='abund.dat',status='unknown')
DO k=1,igt
    write(19,*) ispv(k)
END DO
CLOSE(19)

iptime = ntime
iitime = iitime + iptime
CALL WRITEPOP
CALL WRITEPHYLOGENY


END PROGRAM ibm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE JUMPTABLE(xj,yj)
integer(2) xj(32),yj(32)
xj(1)=0;    yj(1)=1;
xj(2)=1;    yj(2)=0;
xj(3)=0;    yj(3)=-1;
xj(4)=-1;   yj(4)=0;  !first 4 neighbors, dist-max = 1

xj(5)=1;    yj(5)=1;
xj(6)=1;    yj(6)=-1;
xj(7)=-1;   yj(7)=-1;
xj(8)=-1;   yj(8)=1;
xj(9)=0;    yj(9)=2;
xj(10)=2;   yj(10)=0;
xj(11)=0;   yj(11)=-2;
xj(12)=-2;  yj(12)=0;
xj(13)=1;   yj(13)=2;
xj(14)=2;   yj(14)=1;
xj(15)=2;   yj(15)=-1;
xj(16)=1;   yj(16)=-2;
xj(17)=-1;  yj(17)=-2;
xj(18)=-2;  yj(18)=-1;
xj(19)=-2;  yj(19)=1;
xj(20)=-1;  yj(20)=2; !first 20 neighbors, dist-max = 2

xj(21)=0;   yj(21)=3;
xj(22)=3;   yj(22)=0;
xj(23)=0;   yj(23)=-3;
xj(24)=-3;  yj(24)=0;
xj(25)=1;   yj(25)=3;
xj(26)=3;   yj(26)=1;
xj(27)=3;   yj(27)=-1;
xj(28)=1;   yj(28)=-3;
xj(29)=-1;  yj(29)=-3;
xj(30)=-3;  yj(30)=-1;
xj(31)=-3;  yj(31)=1;
xj(32)=-1;  yj(32)=3; !first 32 neighbors, dist-max = 3

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For each individual k, find a kmother and its mate      !
! Mother is usually k, but not always. 
! At the beggining NEIG was called with irad = 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDMATE(k,kmother,kmate,irad)
USE globals
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(OUT) :: kmate
INTEGER, INTENT(INOUT) :: kmother
INTEGER(2), ALLOCATABLE :: neigs(:)
INTEGER(2) ichoose,imother,ineighbork

ALLOCATE (neigs(nc))
CALL RANDOM_NUMBER(aux)

! choose a new mother if aux < Q or if potential mates < P
! the choice may involve increasing irad until new mother has at least 2 potential mates
!
IF(aux < qmat .OR. ineighborg < mnpm) THEN
	CALL RANDOM_NUMBER(aux)
	ichoose = int(aux*ineighbor)+1
	kmother = neigsp(ichoose)  ! choose a new mother
	if(kmother < 1) write(6,*) aux,ichoose
	CALL FINDNEIG(kmother,irad) ! find mother's neighbors
	icc = 0
	DO WHILE (ineighborg < 2 .AND. irad < nrmax)  ! choose a neighbor mother with ineighborg > 1
		IF(icc > 1) THEN
			irad = irad + 1  ! increase search radius
			icc = 0
			CALL FINDNEIG(k,irad)
		END IF
		CALL RANDOM_NUMBER(aux)  !try a different neighbor within the same irad
		ichoose = int(aux*ineighbor)+1
		kmother = neigsp(ichoose)
		! if there are no spatial neighbors kmother=0 and need to increase neighborhood size
		IF(kmother /= 0) THEN
			CALL FINDNEIG(kmother,irad) ! find mother's neighbors
		END IF
		icc = icc + 1
	END DO
END IF


! choose a mate for the mother
!
IF (ineighborg >= 2) THEN
        CALL RANDOM_NUMBER(aux)
	kmate = neig(int(aux*ineighborg)+1)
	IF(kmate == kmother ) then
		write(6,*) 'cannot find a mother'
		STOP
	END IF
END IF

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find spatial and genetic neighbors of kmother !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDNEIG(kmother,irad)
USE globals
INTEGER, INTENT(IN) :: kmother

ineighbor = 0
ineighborg = 0
iibound = 0
neig = 0   ! genetic neighbors
neigsp = 0 ! spatial neighbors 
! the focal individual kmother is not included as a neighbor of itself
ix = x(kmother)
iy = y(kmother)
loop1: DO isite = 1,ndq(irad) 
		ix1 = ix+idqx(irad,isite)
		iy1 = iy+idqy(irad,isite)
		IF (ix1 > nf1) ix1 = ix1 - nf1
		IF (ix1 < 1) ix1 = nf1 + ix1
		IF (iy1 > nf2) iy1 = iy1 - nf2
		IF (iy1 < 1) iy1 = nf2 + iy1
		iibound = iibound + 1 ! number of actual sites in the neighborhood
		loop2: DO iworg = 1,nworg(ix1,iy1)
        			if(worg(ix1,iy1,iworg) /= kmother) then
            				ineighbor = ineighbor + 1
            				neigsp(ineighbor) = worg(ix1,iy1,iworg)
            				dista = 0
            				DO l=1,nb
                				dista = dista + ABS(g(kmother,l)-g(neigsp(ineighbor),l))
                				IF(dista > rg) CYCLE loop2
            				END DO
            				ineighborg = ineighborg + 1
            				neig(ineighborg) = neigsp(ineighbor)
        			end if
		END DO loop2
END DO loop1

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSPECIES
USE globals, ONLY: nf1,nf2,nct,nb,rg,g,igt,ispecies,ispv,x,y,ispidx
INTEGER(4), ALLOCATABLE :: species(:),auxy1(:),auxy2(:)
INTEGER(2) iix,iiy

ALLOCATE (species(nct),auxy1(nct),auxy2(nct))

itot = 0  ! count total population in groups
igt = 0   ! count number of groups
i2 = nct
DO i=1,i2  ! initialize aux2 -- contains all individuals not yet classified
	auxy2(i) = i
END DO 
	
DO WHILE (itot < nct)
	icr = auxy2(1)   !take first individual and find its species
	isp = 0
	ispold = 1
	i1 = 0
	auxy1 = 0
	loop1: DO i=1,i2
			ii = auxy2(i)
			dista = 0
			DO l=1,nb
				IF(g(icr,l) /= g(ii,l)) dista = dista + 1  
				IF(dista > rg) THEN
					i1 = i1 + 1
					auxy1(i1) = ii      !put creatures with dist > rg into aux1
					CYCLE loop1
				END IF
			END DO
			isp = isp + 1
			species(isp) = ii   !collect individuals with dist <= rg from icr
	END DO loop1

	!check if individuals in aux1 have to be included; put the rest in aux2
	itest = 1
	DO WHILE(itest /= 0)
		i2 = 0
		auxy2 = 0
		itest = 0
		isp0 = isp
		IF(i1 /= 0) THEN
			loop2: DO i=1,i1
					DO ji=ispold+1,isp0  
						dista = 0
						DO l=1,nb
							IF(g(auxy1(i),l) /= g(species(ji),l)) dista = dista + 1  !HERE
							IF(dista > rg) EXIT
						END DO
						IF(dista <= rg) THEN
							isp = isp + 1 
							species(isp) = auxy1(i)   ! colect the aux1 individual
							itest = 1                 ! indicates that the process has to be repeated
							CYCLE loop2
						END IF
					END DO
					i2 = i2 + 1
					auxy2(i2) = auxy1(i)  ! put individual in aux2
			END DO loop2
		END IF
		auxy1 = auxy2   ! aux1 contains the creatures not in the species
		i1 = i2
		ispold = isp0
	END DO

	itot = itot + isp    !total number of individuals classified into species
	igt = igt + 1        !number of species

	! save species info
	DO i=1,isp
		ispecies(igt,i) = species(i)
        	ispidx(species(i)) = igt
	END DO
	ispv(igt) = isp          ! number of individuals in species

END DO

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE
