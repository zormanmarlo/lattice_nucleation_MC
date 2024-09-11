MODULE bonded_list_mod
IMPLICIT NONE
PUBLIC :: bonded_list_type
TYPE bonded_list_type
   INTEGER :: lat_vec ( 3 )
END TYPE bonded_list_type
END MODULE bonded_list_mod

PROGRAM RosenBluthMC
USE bonded_list_mod, ONLY : bonded_list_type
! Iteration variables
IMPLICIT NONE
include 'mpif.h'

INTERFACE 
  SUBROUTINE get_bonded_site ( bl_type, latticeXYZ, latticeIDX, jparticle, Nl, tmpPart, rn, XYZ )
       USE bonded_list_mod, ONLY : bonded_list_type
       INTEGER, PARAMETER :: dp =8
       TYPE ( bonded_list_type ), INTENT ( IN ) :: bl_type (:)
       INTEGER, INTENT ( IN ) :: latticeXYZ (:,:)
       INTEGER, INTENT ( IN ) :: latticeIDX (:,:,:)
       INTEGER, INTENT ( IN ) :: jparticle, Nl
       INTEGER, INTENT ( OUT ) :: tmpPart
       REAL ( dp ), INTENT ( IN ) :: rn
       INTEGER, OPTIONAL, INTENT (OUT) :: XYZ (3)
   END SUBROUTINE get_bonded_site 
END INTERFACE

INTERFACE 
   SUBROUTINE swap ( latticeXYZ, latticeIDX, tmpXYZ, tmpInitXYZ, tmpPart, tmpPart2, lswap )
     INTEGER, INTENT ( INOUT ) :: latticeXYZ (:,:)
     INTEGER, INTENT ( INOUT ) :: latticeIDX (:,:,:)
     INTEGER, INTENT ( INOUT ) ::  tmpInitXYZ (:)
     INTEGER, INTENT ( IN ) ::  tmpXYZ (:), tmpPart, tmpPart2
     LOGICAL, INTENT ( IN ) :: lswap
   END SUBROUTINE swap 
END INTERFACE

INTERFACE 
   SUBROUTINE get_energy ( dE, Part, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
      INTEGER, PARAMETER :: dp = 8
      REAL ( dp ), INTENT ( OUT ) :: dE
      INTEGER, INTENT ( IN ) :: Part, Nl, Ntot
      INTEGER, INTENT ( IN ) :: latticeXYZ (:,:), idxType (:)
      REAL ( dp ), INTENT ( IN ) :: fQ(0:,0:,0:), fu (0:,0:,0:,0:,0:)
      REAL ( dp ), INTENT ( IN ) :: QQ, zz(0:), kT
   END SUBROUTINE 
END INTERFACE

INTERFACE 
   FUNCTION md(i,Nl)
    INTEGER :: md
    INTEGER, INTENT (IN)  :: i,Nl
   END FUNCTION md
END INTERFACE

INTERFACE 
   SUBROUTINE reallocate ( p, lb1_new, ub1_new )
     IMPLICIT NONE
     INTEGER, PARAMETER :: dp = 8
     INTEGER, DIMENSION ( : ), POINTER :: p
     INTEGER, INTENT ( IN ) :: lb1_new, ub1_new
   END SUBROUTINE reallocate 
END INTERFACE

INTERFACE 
   SUBROUTINE count(Nl, Ntot, ix, nm, ixm, nam, bin, namx, max0, bin0, targetIdx)
    INTEGER, PARAMETER :: dp = 8
    INTEGER :: Nl, Ntot, nm, namx, targetIdx
    INTEGER :: max0
    INTEGER :: ix(Nl**3, 3), bin(Nl*Nl*Nl)
    INTEGER :: ixm(Ntot, Ntot), nam(Ntot), bin0(Ntot)
   END SUBROUTINE count
END INTERFACE


INTEGER, PARAMETER :: dp = 8
INTEGER :: Nl ! Lattice size
INTEGER :: x,y,z, irb, nrb, rbindex
INTEGER :: i,j,k, ii, idxTrial
INTEGER :: currentPart, tmpPart,tmpPart2,tmpAnchor,tmpInitPart,tmpSwap
INTEGER :: RBPart, SwapPart, icount
REAL ( dp ) :: start_outin, finish_outin
REAL ( dp ) :: start_wnew, finish_wnew
REAL ( dp ) :: start_wold, finish_wold
LOGICAL :: lavbmc_reject, ltrans_reject
INTEGER :: tmpTargetSize, targetSize, currentSize
REAL ( dp ) :: ener, rand_weight, wnew, wold, running_weight
REAL ( dp ), ALLOCATABLE :: w (:), enew ( : ), eold ( : )
INTEGER, ALLOCATABLE :: part_rb ( : )

! lattice translation array
INTEGER  :: LatTrans ( -3:3, -3:3, -3:3 )
INTEGER, ALLOCATABLE  :: BondedList ( : )
REAL ( dp )  :: Vin, Vout
INTEGER :: Nin, Nout, BondedPart, isite, jParticle
INTEGER :: TrialSite
LOGICAL :: lfound

! RanDOm / temporary variables
REAL ( dp ) :: RAN1
INTEGER :: IDUM
REAL ( dp ) :: tmp,tmp2,rand,rn
LOGICAL :: linout, loutin, ltrans

! overflow/underflow
REAL (dp) :: exp_max_val, exp_min_val
REAL (dp) :: min_val, max_val, exp_test

! Loop / output frequency variables
INTEGER :: itraj,ntraj,nequil
INTEGER :: enWriteFreq,xyzWriteFreq

! Number of Ca, CO3, and total number of particles
INTEGER :: N1,N2,Ntot,tmpNewXYZ
INTEGER :: nClusters, maxSize, tmpMaxSize
INTEGER :: tmpX,tmpY,tmpZ,trans

! Testing variables
INTEGER :: nAccept_trans, nTrial_trans
INTEGER :: nAccept_inout, nTrial_inout
INTEGER :: nAccept_outin, nTrial_outin

! Thermal constants
REAL ( dp ) ::  temp
REAL ( dp ) ::  kT

! fQ is coulombic interaction
! fu is PMF interaction
REAL ( dp ), ALLOCATABLE :: fQ(:,:,:),fu(:,:,:,:,:)

! Lattice and particle arrays
INTEGER, ALLOCATABLE :: lattice(:,:,:),latticeXYZ(:,:),idxType(:),latticeIDX(:,:,:)
INTEGER, ALLOCATABLE :: tmpInitXYZ(:), tmpXYZ(:), transXYZ(:)
INTEGER, ALLOCATABLE :: SwapXYZ(:), RBXYZ(:)
INTEGER, ALLOCATABLE :: bin(:), tmpBin(:)
INTEGER, ALLOCATABLE :: allClusterSizes(:), clusterArray(:,:), targetClusterIdx(:), inSpaces(:)
INTEGER, ALLOCATABLE :: clusterList(:), tmpClusterList(:)

! MPI variables
INTEGER :: ierr,np,nw,me
CHARACTER(len=2) fme

! Energy and biasing variables
INTEGER :: nBiasBins
REAL ( dp ) :: Ewrite,Ebiased,EwriteBiased,dEt,dEbiased,Et
REAL ( dp ), ALLOCATABLE :: Kn(:)
REAL ( dp ) :: E,dE
REAL ( dp ) :: zz(0:2),QQ
REAL ( dp ) :: avbmcE
INTEGER :: iostat

! Lattice initialization variables
LOGICAL :: restart
INTEGER :: irestart
CHARACTER(len=200) :: filename
INTEGER :: xPos, yPos, zPos
REAL ( dp ) :: xPos_tmp, yPos_tmp, zPos_tmp
INTEGER :: nParticles, particleType
CHARACTER(len=2) :: element
TYPE (bonded_list_type), ALLOCATABLE :: bl_type (:)


! Initialize MPI so we can run in parallel
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)
! me = 0
write(fme,'(i2.2)') me

! Zero empty-empty
zz(0)=  0.0_dp

! Read bias file
OPEN(10,FILE='Kn')
READ(10,*) nBiasBins
ALLOCATE(Kn(nBiasBins))
DO ii = 1, SIZE ( Kn )
  READ(10,*) Kn(ii)
END DO
CLOSE(10)

! Read input file
OPEN(10,FILE='flp.in')
READ(10,*) Nl,IDUM,ntraj,temp,nequil,enWriteFreq,xyzWriteFreq
READ(10,*) N1,N2
READ(10,*) zz(1),zz(2),QQ,kT
READ(10,*) targetSize, nrb
READ(10,*) irestart
CLOSE(10)
restart = .FALSE.
IF (irestart.eq.1) restart = .TRUE.
Ntot=N1+N2
WRITE ( *, * ) 'NTOT ', Ntot
WRITE ( *, * ) 'Number of RB trials ', nrb
ALLOCATE ( w (nrb ) )
ALLOCATE ( part_rb (nrb ) )
ALLOCATE ( enew (nrb ) )
ALLOCATE ( eold (nrb ) )
! Initialize interaction and position array
ALLOCATE(fQ(0:Nl/2,0:Nl/2,0:Nl/2))
ALLOCATE(fu(0:2,0:2,0:Nl/2,0:Nl/2,0:Nl/2))
ALLOCATE(lattice(Nl,Nl,Nl))
ALLOCATE(latticeXYZ(Nl**3,3))
ALLOCATE(idxType(Nl**3))
ALLOCATE(latticeIDX(Nl,Nl,Nl))

ALLOCATE(tmpXYZ(3))
ALLOCATE(RBXYZ(3))
ALLOCATE(SwapXYZ(3))
ALLOCATE(tmpInitXYZ(3))
ALLOCATE(transXYZ(3))

ALLOCATE(clusterArray(Ntot,Ntot))
ALLOCATE(allClusterSizes(Ntot))
ALLOCATE(targetClusterIdx(Ntot))
ALLOCATE(bin(Nl*Nl*Nl), tmpBin(Nl*Nl*Nl))
ALLOCATE(inSpaces(Nl**3))
ALLOCATE(clusterList(Ntot), tmpClusterList(Ntot))

! Read potentials file
i=0
OPEN(10,FILE='f1')
DO
  READ(10,*,iostat=iostat) x,y,z,fQ(x,y,z),&
      & fu(1,1,x,y,z),&
      & fu(1,2,x,y,z),&
      & fu(2,2,x,y,z)
  if(iostat.lt.0) exit
  fu(0,0,x,y,z)=0.0_dp
  fu(0,1,x,y,z)=0.0_dp
  fu(0,2,x,y,z)=0.0_dp
  fu(1,0,x,y,z)=0.0_dp
  fu(2,0,x,y,z)=0.0_dp
  fu(2,1,x,y,z)=fu(1,2,x,y,z)
  i=i+1
ENDDO
CLOSE(10)

! Generate ID for MPI thREAD and OPEN output files
IDUM=-abs(IDUM+me)
OPEN(11,FILE='E-'//fme//'.log')
OPEN(12,FILE='target_cluster-'//fme//'.out')
OPEN(13,FILE='traj-'//fme//'.xyz')
OPEN(14,FILE='stats-'//fme//'.log')
OPEN(15,FILE='clusters-'//fme//'.out')

! Initialize random variable and cluster lists
tmp=RAN1(IDUM)

! define overflow/underflow bounds
exp_max_val=0.9_dp*LOG (HUGE(0.0_dp))
exp_min_val=0.9_dp*LOG (TINY(0.0_dp))
max_val=EXP(exp_max_val)
min_val=0.0_dp
WRITE ( *, * ) 'MAX VAL ', exp_max_val, max_val
WRITE ( *, * ) 'MIN VAL ', exp_min_val, min_val

! Initialize array
lattice(:,:,:)=0

! Construct lattice translations that are within cuttoff of central site
LatTrans (:,:,:) = 0
DO i = -3, 3
  DO j = -3, 3
    DO k = -3,3
      IF ( i == 0 .AND. j==0 .AND. k==0 ) CYCLE
      IF ( REAL (i**2, dp)+ REAL (j**2, dp) + REAL (k**2, dp) <= 2.0_dp) &
        LatTrans ( i, j, k ) = 1
    END DO
  END DO
END DO
Vin = REAL ( SUM ( LatTrans ), dp )
Vout = REAL (Nl**3, dp)  - Vin
WRITE ( *,  * ) 'Number of Bonded Sites = ', SUM (LatTrans) 
ALLOCATE (BondedList ( SUM (LatTrans ) ) )
ALLOCATE ( bl_type ( SUM(LatTrans) ) )
icount = 0
DO i=-3,3
  DO j=-3,3
    DO k=-3,3
      IF (LatTrans (i,j,k)==1 ) THEN
        icount=icount + 1
        bl_type (icount)%lat_vec (:) = (/i,j,k/)
      END IF
    END DO
  END DO
END DO
WRITE ( *, * ) 'ICOUNT ', icount, 'AVAL Lattice Translation ', SUM(LatTrans) 

! Initialize particle array
! Either by reading in input file or generating random positions
IF (restart) THEN
  filename = "input-" // trim(fme) // ".xyz"
  OPEN(unit=10, file=TRIM(filename), status='old', action='READ')
  READ(10, *) nParticles
  READ(10, *)
  ! Loop over each line to read and set the particle data
  DO i = 1, nParticles
      READ(10, *) element, xPos_tmp, yPos_tmp, zPos_tmp
      IF (TRIM(element) == 'H') THEN
          particleType = 1
      ELSE IF (TRIM(element) == 'O') THEN
          particleType = 2
      ELSE
          WRITE(*, *) 'Unknown element type in file'
          STOP
      END IF
      xPos = INT(xPos_tmp)
      yPos = INT(yPos_tmp)
      zPos = INT(zPos_tmp)
      latticeXYZ(i, :) = [xPos, yPos, zPos]
      idxType(i) = particleType
      lattice(xPos, yPos, zPos) = i
      latticeIDX(xPos, yPos, zPos) = i
  END DO
  CLOSE(10)
ELSE
  DO i=1,N1
    DO
      xPos=INT(REAL(Nl,dp)*RAN1(IDUM))+1
      yPos=INT(REAL(Nl,dp)*RAN1(IDUM))+1
      zPos=INT(REAL(Nl,dp)*RAN1(IDUM))+1
      IF (lattice(xPos,yPos,zPos)==0) EXIT
    ENDDO
    idxType(i)= 1
    latticeXYZ(i,1)=xPos
    latticeXYZ(i,2)=yPos
    latticeXYZ(i,3)=zPos
    lattice(xPos,yPos,zPos)=i
    latticeIDX(xPos,yPos,zPos) = i
  ENDDO
  DO i=N1+1,Ntot
    DO
      xPos=int(REAL(Nl,dp)*RAN1(IDUM))+1
      yPos=int(REAL(Nl,dp)*RAN1(IDUM))+1
      zPos=int(REAL(Nl,dp)*RAN1(IDUM))+1
      IF(lattice(xPos,yPos,zPos)==0) EXIT
    ENDDO
    idxType(i)= 2
    latticeXYZ(i,1)=xPos
    latticeXYZ(i,2)=yPos
    latticeXYZ(i,3)=zPos
    lattice(xPos,yPos,zPos)=i
    latticeIDX(xPos,yPos,zPos) = i
  ENDDO
ENDIF
i=Ntot
DO xPos=1,Nl
  DO yPos=1,Nl
    DO zPos=1,Nl
      IF(lattice(xPos,yPos,zPos)==0)THEN
        i=i+1
        idxType(i)=  0
        latticeXYZ(i,1)=xPos
        latticeXYZ(i,2)=yPos
        latticeXYZ(i,3)=zPos
        lattice(xPos,yPos,zPos)=i
        latticeIDX(xPos,yPos,zPos) = i
      ENDIF
    ENDDO
  ENDDO
ENDDO

! Calculate energy of initial configuration
E=0.0_dp
Ebiased=0.0_dp
DO i=1,Ntot
  DO j=1,Ntot
    IF (i/=j) THEN
      x=md(latticeXYZ(i,1)-latticeXYZ(j,1),Nl)
      y=md(latticeXYZ(i,2)-latticeXYZ(j,2),Nl)
      z=md(latticeXYZ(i,3)-latticeXYZ(j,3),Nl)
      E=E+0.5_dp*fQ(x,y,z)*zz(idxType(i))*zz(idxType(j))*QQ
      E=E+0.5_dp*kT*fu(idxType(i),idxType(j),x,y,z)
    ENDIF
  ENDDO
ENDDO

! Initalize cluster settings
targetClusterIdx(:) = 0
targetClusterIdx(1) = 1
allClusterSizes(:) = 1
clusterList(:) = 0
clusterList(1) = 1
tmpClusterList = clusterList
!targetSize = 1
maxSize = 1
tmpMaxSize = 1
CALL count(Nl, Ntot, latticeXYZ, nClusters, clusterArray, allClusterSizes, &
          bin, maxSize, currentSize, clusterList, 1)

     WRITE ( *, * ) 'Current Size of Tagged Cluster ', currentSize
     nAccept_trans = 0
     nAccept_inout = 0
     nAccept_outin = 0
     nTrial_outin = 0
     nTrial_inout = 0
     nTrial_trans = 0
! Iterate over nequil+ntraj steps
     
    ! IF (ASSOCIATED ( BondedList)) NULLIFY ( BondedList )
    ! IF (.NOT. ASSOCIATED (BondedList ) ) ALLOCATE ( BondedList ( 1:1 ) )
     DO itraj= -nequil,ntraj
        linout = .FALSE.
        loutin = .FALSE.
        ltrans = .FALSE.
! Randomly select either translation or AVBMC move
       tmp=RAN1(IDUM)
       IF (tmp  <= 0.25_dp) THEN  ! Perform AVBMC-2 Move
         lavbmc_reject = .FALSE.
!Randomly select particle j
         tmp=RAN1(IDUM)
         jParticle=INT(REAL(Ntot,dp)*tmp)+1
! create bonded list of jParticle
         Nin = 0
         DO i = -3, 3
           DO j = -3, 3
             DO k = -3, 3
               IF (LatTrans (i,j,k) == 1 ) THEN
                 transXYZ ( : ) = (/i,j,k/)
                 DO ii=1,3
                   IF ((latticeXYZ(jParticle,ii)+transXYZ(ii) )>Nl) THEN
                     tmpXYZ(ii) = transXYZ(ii)
                   ELSE IF (latticeXYZ(jParticle,ii)+transXYZ(ii)<1) THEN
                     tmpXYZ(ii) = Nl + transXYZ(ii)
                   ELSE
                     tmpXYZ(ii) = latticeXYZ(jParticle,ii)+transXYZ(ii)
                   ENDIF
                 ENDDO
                 BondedPart = latticeIDX(tmpXYZ(1),tmpXYZ(2),tmpXYZ(3))
                 IF ( idxType(BondedPart) /=  0 ) THEN 
                   Nin = Nin + 1
    !               IF ( Nin > SIZE ( BondedList ) ) CALL reallocate (BondedList, 1, Nin )
                   BondedList ( Nin ) = BondedPart
                 END IF
               END IF
             END DO
           END DO
         END DO
         Nout = Ntot - Nin
         tmp2=RAN1(IDUM)
         IF ( Nin == 0 ) tmp2=0.0_dp ! necessitate an out-->in move
         IF ( tmp2 <= 0.5_dp ) THEN !  out-->in
           loutin = .TRUE.
           !CALL CPU_TIME ( start_outin )
           nTrial_outin = nTrial_outin + 1
           ! find a particle in the out region of j 
           DO
             tmpPart2=INT(REAL(Ntot,dp)*RAN1(IDUM))+1
             tmpXYZ(:)=latticeXYZ(tmpPart2,:)
             x=md(tmpXYZ(1)-latticeXYZ(jParticle,1),Nl)
             y=md(tmpXYZ(2)-latticeXYZ(jParticle,2),Nl)
             z=md(tmpXYZ(3)-latticeXYZ(jParticle,3),Nl)
             IF ( REAL(x**2,dp) + REAL(y**2,dp) + REAL(z**2,dp) > 2.0_dp ) EXIT 
           ENDDO
           SwapXYZ = tmpXYZ  ! save the coordinates of the swaping particle
           SwapPart=tmpPart2 ! save label of SwapPart
           wnew = 0.0_dp
           w(:) = 0.0_dp
           !CALL CPU_TIME ( start_wnew )
           DO irb = 1, nrb  ! pick n-sites and compute their rosenbluth weights
           ! pick a UNIQUE random site in the bonded region of jParticle
             DO
               lfound = .FALSE.
               rn = RAN1(IDUM)
               CALL get_bonded_site ( bl_type, latticeXYZ, latticeIDX, jparticle, Nl, RBPart, rn )
               part_rb ( irb ) = RBPart
               ! Screen to see if has been chosen
               DO ii = 1, irb-1
                 IF ( part_rb ( ii ) == RBpart ) lfound = .TRUE.
                 IF ( lfound ) EXIT
               END DO
               IF ( .NOT. lfound ) EXIT
             END DO
             IF ( idxType(RBPart) /=  0 ) THEN 
               w ( irb ) = min_val
             ELSE ! compute rb weight
               dE=0.0_dp
               ! swap SwapPart with RBPart
               CALL swap ( latticeXYZ, latticeIDX, SwapXYZ, RBXYZ, RBPart, SwapPart, .TRUE. )
               ! compute energy of SwapPart in the bonded region
               CALL get_energy ( ener, SwapPart, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
               dE=dE + ener
               !Enew (irb) = E + dE
               Enew (irb) = dE
               exp_test = -Enew (irb)/kT 
               IF ( exp_test > exp_max_val ) THEN
                 w ( irb ) = max_val
               ELSEIF ( exp_test < exp_min_val ) THEN
                 w ( irb ) = min_val
               ELSE
                 w ( irb ) = exp ( exp_test )
               END IF
               ! swap  RBPart with SwapPart to original configuration
               CALL swap ( latticeXYZ, latticeIDX, SwapXYZ, RBXYZ, RBPart, SwapPart, .FALSE. )
             END IF
             wnew = wnew + w (irb)
           END DO !rb loop
           !CALL CPU_TIME ( finish_wnew )
        !   WRITE ( *, * ) 'TIME wnew' , finish_wnew-start_wnew
  ! draw a configuration from rb trials
           rand_weight =  wnew*RAN1(IDUM) 
           running_weight = 0.0_dp
           DO irb = 1, nrb
             running_weight= running_weight+ w ( irb )
             IF (running_weight > rand_weight) THEN 
                RBPart = part_rb ( irb )
                rbindex = irb
                EXIT
             END IF
           END DO 
           w(:)=0.0_dp
           wold=0.0_dp
           dE=0.0_dp
           !compute Energy of SwapPart in original state
           CALL get_energy ( ener, SwapPart,  latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
           dE=dE + ener
           !Eold ( 1 ) = E + dE
           Eold ( 1 ) = dE
           exp_test = -Eold ( 1 )/kT 
           IF ( exp_test > exp_max_val ) THEN
              w ( 1 ) = max_val
           ELSEIF ( exp_test < exp_min_val ) THEN
              w ( 1 ) = min_val
           ELSE
              w ( 1 ) = exp ( exp_test )
           END IF
           wold = wold + w ( 1 )
           !CALL CPU_TIME ( start_wold )
           DO irb = 2, nrb  ! pick n sites and compute their rosenbluth weights
           ! find space in out region to swap
             DO
               tmpXYZ(1)=INT(REAL(Nl,dp)*RAN1(IDUM))+1
               tmpXYZ(2)=INT(REAL(Nl,dp)*RAN1(IDUM))+1
               tmpXYZ(3)=INT(REAL(Nl,dp)*RAN1(IDUM))+1
               x=md(tmpXYZ(1)-latticeXYZ(jParticle,1),Nl)
               y=md(tmpXYZ(2)-latticeXYZ(jParticle,2),Nl)
               z=md(tmpXYZ(3)-latticeXYZ(jParticle,3),Nl)
               tmpPart2=latticeIDX(tmpXYZ(1),tmpXYZ(2),tmpXYZ(3))
               IF ( idxTYPE(tmpPart2) /= 0 ) CYCLE
               IF ( REAL(x**2,dp) + REAL(y**2,dp) + REAL(z**2,dp) > 2.0_dp ) EXIT 
             ENDDO
             !CALL CPU_TIME ( finish_wold )
           !  WRITE ( *, * ) 'TIME wold' , finish_wold-start_wold
             dE=0.0_dp
             CALL swap ( latticeXYZ, latticeIDX, SwapXYZ, tmpInitXYZ, tmpPart2, SwapPart, .TRUE. )
             CALL get_energy ( ener, SwapPart,  latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
             dE=dE + ener
             CALL swap ( latticeXYZ, latticeIDX, SwapXYZ, tmpInitXYZ, tmpPart2, SwapPart, .FALSE. )
             Eold(irb)=dE
             exp_test = -Eold (irb)/kT 
             IF ( exp_test > exp_max_val ) THEN
               w ( irb ) = max_val
             ELSEIF ( exp_test < exp_min_val ) THEN
               w ( irb ) = min_val
             ELSE
               w ( irb ) = exp ( exp_test )
             END IF
             wold = wold + w ( irb )
           END DO
           !CALL CPU_TIME ( finish_outin )
          ! WRITE ( *, * ) 'TIME out-->in' , finish_outin-start_outin
         ELSE  ! in-->out
           linout = .TRUE.
           nTrial_inout = nTrial_inout + 1
           ! find particle in bonded list
           tmpPart = BondedList ( INT(REAL(Nin,dp)*RAN1(IDUM))+1 )
           ! find space in out region to swap
           DO
             tmpXYZ(1)=INT(REAL(Nl,dp)*RAN1(IDUM))+1
             tmpXYZ(2)=INT(REAL(Nl,dp)*RAN1(IDUM))+1
             tmpXYZ(3)=INT(REAL(Nl,dp)*RAN1(IDUM))+1
             x=md(tmpXYZ(1)-latticeXYZ(jParticle,1),Nl)
             y=md(tmpXYZ(2)-latticeXYZ(jParticle,2),Nl)
             z=md(tmpXYZ(3)-latticeXYZ(jParticle,3),Nl)
             IF ( REAL(x**2,dp) + REAL(y**2,dp) + REAL(z**2,dp) > 2.0_dp ) EXIT 
           ENDDO
           tmpPart2=latticeIDX(tmpXYZ(1),tmpXYZ(2),tmpXYZ(3))
           IF ( idxTYPE(tmpPart2) /= 0 ) lavbmc_reject=.TRUE.
         ENDIF
  ! reject move early if possible
         IF (lavbmc_reject ) GOTO 500

         IF ( linout ) THEN ! compute dE for in-->out and swap config
  ! Zero out energy change variables
           dE=0.0_dp
           dEbiased=0.0_dp
  ! calculate energy between target and all other particles
           CALL get_energy ( ener, tmpPart, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
           dE=dE-ener
           CALL get_energy ( ener, tmpPart2, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
           dE=dE-ener
           CALL swap ( latticeXYZ, latticeIDX, tmpXYZ, tmpInitXYZ, tmpPart, tmpPart2, .TRUE. )
  ! Recalculate energy between target and all other particles
           CALL get_energy ( ener, tmpPart, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
           dE=dE+ener
           CALL get_energy ( ener, tmpPart2, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
           dE=dE+ener
         ELSEIF  (loutin) THEN ! swap config
           CALL swap ( latticeXYZ, latticeIDX, SwapXYZ, RBXYZ, RBPart, SwapPart, .TRUE. )
         ELSE
           WRITE (*, *) "bad"
           STOP
         END IF
         CALL count(Nl, Ntot, latticeXYZ, nClusters, clusterArray, allClusterSizes, &
          tmpBin, tmpMaxSize, tmpTargetSize, tmpClusterList, 1)
         dEbiased = 0.0_dp
         dEbiased = REAL(Kn(tmpTargetSize) - Kn(currentSize), dp)

  ! AVBMC-2 acceptance criteria
         IF (loutin) THEN
    ! out-->in acceptance criteria 
           avbmcE = exp(-dEbiased/kT)*wnew/wold*(Vin*REAL(Nout,dp))/(Vout*REAL(Nin+1,dp))
         ELSEIF (linout) THEN
    ! in-->out acceptance criteria
           dEt=dE+dEbiased
           avbmcE = exp(-dEt/kT)*(Vout*REAL(Nin,dp))/(Vin*REAL(Nout+1,dp))
         ENDIF
         IF ( avbmcE >= 1.0_dp ) THEN
           rand = 0.0_dp
         ELSE 
           rand = RAN1(IDUM)
         END IF
         IF (rand < avbmcE ) THEN ! Accept
           IF (loutin) nAccept_outin = nAccept_outin + 1
           IF (loutin) E=E+Enew(rbindex)-Eold(1)
           IF (linout) nAccept_inout = nAccept_inout + 1
           IF (linout) E=E+dE
           bin = tmpBin
           maxSize = tmpMaxSize
           currentSize = tmpTargetSize
           clusterList = tmpClusterList
         ELSE ! swap back
           IF (loutin) CALL swap ( latticeXYZ, latticeIDX, SwapXYZ, RBXYZ, RBPart, SwapPart, .FALSE. )
           IF (linout) CALL swap ( latticeXYZ, latticeIDX, tmpXYZ, tmpInitXYZ, tmpPart, tmpPart2, .FALSE. )
         ENDIF
       ELSE !Translations
         nTrial_trans = nTrial_trans + 1
         ltrans_reject = .FALSE.
         ltrans=.TRUE.
  ! Pick particle to move and get translated positions
         tmpPart=INT(REAL(Ntot,dp)*RAN1(IDUM))+1
  ! Pick a random site in the bonded region of tmpPart
         rn = RAN1 (IDUM )
         CALL get_bonded_site ( bl_type, latticeXYZ, latticeIDX, tmpPart, Nl, tmpPart2, rn, XYZ=tmpXYZ )
         IF (idxType(tmpPart2)/= 0 ) ltrans_reject=.TRUE. ! site occupied
         IF (ltrans_reject) GOTO 500
  ! Zero out energy change variables
         dE=0.0_dp

  ! Calculate energy between target and all other particles
         CALL get_energy ( ener, tmpPart, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
         dE=dE-ener
         CALL get_energy ( ener, tmpPart2, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
         dE=dE-ener
         CALL swap ( latticeXYZ, latticeIDX, tmpXYZ, tmpInitXYZ, tmpPart, tmpPart2, .TRUE. )
  ! Recalculate energy between target and all other particles
         CALL get_energy ( ener, tmpPart, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
         dE=dE+ener
         CALL get_energy ( ener, tmpPart2, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
         dE=dE+ener
  ! Update clusters with candidate move
        CALL count(Nl, Ntot, latticeXYZ, nClusters, clusterArray, allClusterSizes, &
                tmpBin, tmpMaxSize, tmpTargetSize, tmpClusterList, 1)
        dEbiased=0.0_dp
        dEbiased = REAL(Kn(tmpTargetSize) - Kn(currentSize),dp)

        dEt=dE+dEbiased

  ! Metropolis acceptance criteria
         tmp = RAN1 (IDUM)
         IF( (dEt.lt.0.0_dp).or.(exp(-dEt/temp).gt.tmp))THEN
           E=E+dE
           Ebiased = Ebiased + dEbiased
           nAccept_trans = nAccept_trans + 1
           bin = tmpBin
           maxSize = tmpMaxSize
           currentSize = tmpTargetSize
           clusterList = tmpClusterList
         ELSE
    ! Swap back
           CALL swap ( latticeXYZ, latticeIDX, tmpXYZ, tmpInitXYZ, tmpPart, tmpPart2, .FALSE. )
         ENDIF
       ENDIF

500   CONTINUE
! Write out energy and cluster size
IF (modulo(itraj-1,enWriteFreq)==0)THEN
  write(11,'(3e20.9)') E, Ebiased, E+Ebiased
  CALL flush(11)
  CALL count(Nl, Ntot, latticeXYZ, nClusters, clusterArray, allClusterSizes, &
             bin, maxsize, tmpTargetSize, ClusterList, 1)
  write(12, * ) tmpTargetSize, targetSize
  CALL flush(12)
  write(14, *) nTrial_inout,nAccept_inout, &
               nTrial_outin,nAccept_outin, &
               nTrial_trans,nAccept_trans
  CALL flush(14)
  write(15,'(1000i5)') maxSize,bin(1:maxSize)
  CALL flush(15)
ENDIF

! Write out trajectory
if(modulo(itraj-1,xyzWriteFreq).eq.0)then
  write(13,'(i5)') Ntot
  write(13,'(a,f10.5)') ' temp= ',temp
  DO i=1,N1
    write(13,'(a,3f20.10)') 'H ',(REAL(latticeXYZ(i,k),dp),k=1,3)
  ENDDO
  DO i=N1+1,Ntot
    write(13,'(a,3f20.10)') 'O ',(REAL(latticeXYZ(i,k),dp),k=1,3)
  ENDDO
  CALL flush(12)
  ! Check if any rows in latticeXYZ are the same
  DO i = 1, Ntot-1
    DO j = i+1, Ntot
      if (all(latticeXYZ(i, :) == latticeXYZ(j, :))) then
        print *, "Particles ", i, " and ", j, " in frame ", itraj, " are overlapping"
      endif
    ENDDO
  ENDDO
endif

! End of trajectory DO loop
ENDDO
!print *, nAccept

! Write final configuration to XYZ file
write(13,'(i5)') Ntot
write(13,'(a,f10.5)') ' temp= ',temp
DO i=1,N1
  write(13,'(a,3f20.10)') 'H ',(REAL(latticeXYZ(i,k),dp),k=1,3)
ENDDO
DO i=N1+1,Ntot
  write(13,'(a,3f20.10)') 'O ',(REAL(latticeXYZ(i,k),dp),k=1,3)
ENDDO
CALL flush(13)

! Close output files and finalize MPI
CLOSE(11)
CLOSE(12)
CLOSE(13)
CLOSE(14)
CLOSE(15)
CALL MPI_FINALIZE(ierr)
END PROGRAM RosenBluthMC

REAL ( 8 ) FUNCTION RAN1(IDUM)
save
INTEGER :: idum,m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3,iff
INTEGER :: partXYZ1,partXYZ2,partXYZ3,j
REAL ( 8 ) ::   r,rm1,rm2
DIMENSION R(97)
PARAMETER (M1=139968,IA1=3877,IC1=29573,RM1=1.0d0/M1)
PARAMETER (M2=214326,IA2=3613,IC2=45289,RM2=1.0d0/M2)
PARAMETER (M3=714025,IA3=1366,IC3=150889)
DATA IFF /0/
IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
  IFF=1
  partXYZ1=MOD(IC1-IDUM,M1)
  partXYZ1=MOD(IA1*partXYZ1+IC1,M1)
  partXYZ2=MOD(partXYZ1,M2)
  partXYZ1=MOD(IA1*partXYZ1+IC1,M1)
  partXYZ3=MOD(partXYZ1,M3)
  DO J=1,97
    partXYZ1=MOD(IA1*partXYZ1+IC1,M1)
    partXYZ2=MOD(IA2*partXYZ2+IC2,M2)
    R(J)=(FLOAT(partXYZ1)+FLOAT(partXYZ2)*RM2)*RM1
  ENDDO
  IDUM=1
ENDIF
partXYZ1=MOD(IA1*partXYZ1+IC1,M1)
partXYZ2=MOD(IA2*partXYZ2+IC2,M2)
partXYZ3=MOD(IA3*partXYZ3+IC3,M3)
J=1+(97*partXYZ3)/M3
IF(J.GT.97.OR.J.LT.1) STOP
RAN1=R(J)
R(J)=(dble(partXYZ1)+dble(partXYZ2)*RM2)*RM1
RETURN
END

FUNCTION md(i,Nl)
INTEGER :: md
INTEGER, INTENT (IN) :: i,Nl
md=MOD(i,Nl)
IF(md<0)THEN
  md=md+Nl
ENDIF
IF(md>Nl/2)THEN
  md=Nl-md
ENDIF
RETURN
END

FUNCTION md_full(i,Nl)
  INTEGER :: md_full, i, Nl
  md_full = MOD(i, Nl)
  IF (md_full <= 0) THEN
    md_full = md_full + Nl
  ENDIF
  RETURN
END

SUBROUTINE count(Nl, Ntot, ix, nm, ixm, nam, bin, namx, max0, bin0, targetIdx)
    ! Declare variables
    INTEGER, PARAMETER :: dp = 8
    INTEGER :: Nl, Ntot, nm, namx, targetIdx
    INTEGER :: max0
    INTEGER :: ix(Nl**3, 3), bin(Nl*Nl*Nl)
    INTEGER :: ixm(Ntot, Ntot), nam(Ntot), bin0(Ntot)
    INTEGER :: iam, im
    INTEGER :: i, j, k
    INTEGER :: i1, i2, i3
    INTEGER, ALLOCATABLE :: nf(:), ac(:,:)

    ! Allocate memory for helper arrays
    ALLOCATE(nf(Ntot), ac(Ntot, Ntot))

    ! Initialize adjacency matrix to zero
    ac(:,:) = 0
    ! Loop over all pairs of points
    DO i = 1, Ntot
        DO j = 1, Ntot
            if (i .ne. j) then
                ! Calculate periodic boundary conditions
                i1 = md(ix(i, 1) - ix(j, 1), Nl)
                i2 = md(ix(i, 2) - ix(j, 2), Nl)
                i3 = md(ix(i, 3) - ix(j, 3), Nl)
                ! Check if points are within distance 2 (since sqrt(4) = 2)
                if ((i1*i1 + i2*i2 + i3*i3) <= 2.0_dp) then
                    ! Set adjacency matrix to 1 for connected points
                    ac(i, j) = 1
                    ac(j, i) = 1
                endif
            endif
        ENDDO
    ENDDO

    ! Initialize the first index of the labeling function
    nf(1) = 1
    ! Labeling loop to identify connected components
    DO i = 2, Ntot
        nf(i) = i
        DO j = 1, i - 1
            nf(j) = nf(nf(j))
            if (ac(i, j) == 1) nf(nf(nf(j))) = i
        ENDDO
    ENDDO
    ! Flatten the labeling array to root labels
    DO i = 1, Ntot
        nf(i) = nf(nf(i))
    ENDDO

    ! Initialize cluster name array and index matrix
    nam(:)   = 0
    ixm(:,:) = 0
    nm = 0
    ! Loop to build clusters based on labels
    DO i = 1, Ntot
        DO im = 1, nm
            DO iam = 1, nam(im)
                j = ixm(iam, im)
                if (nf(i) == nf(j)) then
                    nam(im) = nam(im) + 1
                    ixm(nam(im), im) = i
                    goto 102
                endif
            ENDDO
        ENDDO
        nm = nm + 1
        nam(nm) = 1
        ixm(1, nm) = i
102 CONTINUE
    ENDDO

    ! Initialize bin counts and find maximum cluster size
    bin(:) = 0
    bin0(:) = 0
    namx = 0
    max0 = 0
    DO im = 1, nm
      if (nam(im) > namx) namx = nam(im)
      bin(nam(im)) = bin(nam(im)) + 1
      ! Check if the ion at index 1 is in the current cluster and update max0
      DO iam = 1, nam(im)
        IF (ixm(iam, im) == targetIdx) THEN
          max0 = nam(im)
          DO j = 1, nam(im)
            bin0(ixm(j, im)) = 1
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    ! DeALLOCATE dynamic arrays
    DEALLOCATE(nf, ac)

    return
END SUBROUTINE
!-------
        SUBROUTINE get_bonded_site ( bl_type, latticeXYZ, latticeIDX, jparticle, Nl, tmpPart, rn, XYZ )
             USE bonded_list_mod, ONLY : bonded_list_type
             INTEGER, PARAMETER :: dp =8
             TYPE ( bonded_list_type ), INTENT ( IN ) :: bl_type (:)
             INTEGER, INTENT ( IN ) :: latticeXYZ (:,:)
             INTEGER, INTENT ( IN ) :: latticeIDX (:,:,:)
             INTEGER, INTENT ( IN ) :: jparticle, Nl
             INTEGER, INTENT ( OUT ) :: tmpPart
             REAL ( dp ), INTENT ( IN ) :: rn
             INTEGER, OPTIONAL, INTENT (OUT) :: XYZ (3)
             INTEGER:: TrialSite, isite, i,j,k
             INTEGER:: transXYZ(3), myXYZ ( 3 )
  ! find the random site in bonded region of jParticle
             TrialSite=INT(SIZE(bl_type)*rn)+1
  ! find the lattice coordinates of random site in bonded region of jParticle
             transXYZ (:) = bl_type (TrialSite)%lat_vec (:)
             DO k=1,3
               IF ((latticeXYZ(jParticle,k)+transXYZ(k) )>Nl) THEN
                 myXYZ(k) = transXYZ(k)
               ELSE IF ((latticeXYZ(jParticle,k)+transXYZ(k))<1) THEN
                 myXYZ(k) = Nl + transXYZ(k)
               ELSE
                 myXYZ(k) = latticeXYZ(jParticle,k)+transXYZ(k)
               ENDIF
             ENDDO
             tmpPart = latticeIDX(myXYZ(1),myXYZ(2),myXYZ(3))
             IF (PRESENT (XYZ) ) XYZ = myXYZ
         END SUBROUTINE get_bonded_site
!-------------
       SUBROUTINE get_energy ( dE, Part, latticeXYZ, idxType, Nl, Ntot, fQ, fu, zz, QQ, kT )
         INTEGER, PARAMETER :: dp = 8
         REAL ( dp ), INTENT ( OUT ) :: dE
         INTEGER, INTENT ( IN ) :: Part, Nl, Ntot
         INTEGER, INTENT ( IN ) :: latticeXYZ (:,:), idxType (:)
         REAL ( dp ), INTENT ( IN ) :: fQ(0:,0:,0:), fu (0:,0:,0:,0:,0:)
         REAL ( dp ), INTENT ( IN ) :: QQ, zz(0:), kT
         INTEGER :: i, x, y, z
         dE=0.0_dp
         DO i=1,Ntot
           x=md(latticeXYZ(Part,1)-latticeXYZ(i,1),Nl)
           y=md(latticeXYZ(Part,2)-latticeXYZ(i,2),Nl)
           z=md(latticeXYZ(Part,3)-latticeXYZ(i,3),Nl)
           dE=dE+fQ(x,y,z)*zz(idxType(Part))*zz(idxType(i))*QQ
           dE=dE+kT*fu(idxType(Part),idxType(i),x,y,z)
         ENDDO
       END SUBROUTINE 
!-------
       SUBROUTINE swap ( latticeXYZ, latticeIDX, tmpXYZ, tmpInitXYZ, tmpPart, tmpPart2, lswap )
         INTEGER, INTENT ( INOUT ) :: latticeXYZ (:,:)
         INTEGER, INTENT ( INOUT ) :: latticeIDX (:,:,:)
         INTEGER, INTENT ( INOUT ) ::  tmpInitXYZ (:)
         INTEGER, INTENT ( IN ) ::  tmpXYZ (:), tmpPart, tmpPart2
         LOGICAL, INTENT ( IN ) :: lswap
         INTEGER :: k
         IF ( lswap ) THEN
  !swap: saving current location in case we need to swap back
           DO k=1,3
             tmpInitXYZ(k)=latticeXYZ(tmpPart,k) ! save current xyz of tmpPart
             latticeXYZ(tmpPart,k)=tmpXYZ(k)     ! move xyz of tmpPart to xyz of tmpPart2 
             latticeXYZ(tmpPart2,k)=tmpInitXYZ(k) ! move xyz of tmpPart2 to xyz of tmpPart
           ENDDO
           latticeIDX(tmpXYZ(1),tmpXYZ(2),tmpXYZ(3)) = tmpPart
           latticeIDX(tmpInitXYZ(1),tmpInitXYZ(2),tmpInitXYZ(3)) = tmpPart2

! Swap back
         ELSEIF (.NOT. lswap ) THEN 
           DO k=1,3
             latticeXYZ(tmpPart,k)=tmpInitXYZ(k)
             latticeXYZ(tmpPart2, k) = tmpXYZ(k)
           ENDDO
           latticeIDX(tmpXYZ(1),tmpXYZ(2),tmpXYZ(3)) = tmpPart2
           latticeIDX(tmpInitXYZ(1),tmpInitXYZ(2),tmpInitXYZ(3)) = tmpPart
         ELSE
           WRITE ( *, * ) 'Error'
           STOP
         ENDIF 
       END SUBROUTINE swap
!-------------
       SUBROUTINE reallocate ( p, lb1_new, ub1_new )
         INTEGER, PARAMETER :: dp = 8
         INTEGER, DIMENSION ( : ), POINTER :: p
         INTEGER, INTENT ( IN ) :: lb1_new, ub1_new
         INTEGER, PARAMETER :: zero=0
         INTEGER, ALLOCATABLE, DIMENSION ( : ) :: work
         INTEGER :: istat, lb1, lb1_old, ub1, ub1_old


         IF (ASSOCIATED(p)) THEN
          lb1_old = LBOUND(p,1)
          ub1_old = UBOUND(p,1)
          ALLOCATE (work(lb1_old:ub1_old),STAT=istat)
          IF (istat /= 0) WRITE (*, * ) "ERROR ALLOCATING work"
          work(:) = p(:)
          DEALLOCATE (p,STAT=istat)
          IF (istat /= 0) WRITE (*, * ) "ERROR DEALLOCATING p"
        END IF


        ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
        IF (istat /= 0)  WRITE (*,*) "ERROR in ALLOCATING p"
        p(:) = zero

        IF (ALLOCATED(work)) THEN
          lb1 = MAX(lb1_new,lb1_old)
          ub1 = MIN(ub1_new,ub1_old)
          p(lb1:ub1) = work(lb1:ub1)
          DEALLOCATE (work,STAT=istat)
          IF (istat /= 0) WRITE ( *, * ) "ERROR in DEALLOCATING work"
        END IF
       END SUBROUTINE reallocate
