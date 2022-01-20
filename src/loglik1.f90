
! module mla1
!   implicit none
!   integer,save ::ny,ns,nv,idiag,nvc,nea,ncor &
!        ,maxmes,nobs,nef,ncontr,ntrtot,nalea &
!        ,nySPL,ntotvalSPL,nyORD,ntotvalORD,npmtot &
!        ,nMC,methInteg,nmescur &
!        ,nvarxevt,nbevt,logspecif,idtrunc,nvdepsurv,nrisqtot,nxevt,nasso
!   integer,dimension(:),allocatable,save::typrisq,nz,nprisq,nevtparx,nxcurr
!   double precision,dimension(:),allocatable,save::Y,uniqueY,minY,maxY,rangeY
!   double precision,dimension(:,:),allocatable,save ::X,zi
!   double precision,dimension(:),allocatable,save::Tsurv0,Tsurv,Tsurvint
!   double precision,dimension(:,:),allocatable,save::Tsurv0_st2,Tsurv_st2
!   integer,dimension(:),allocatable,save::Devt,ind_survint
!   integer,dimension(:),allocatable,save::idtdv,idsurv
!   integer,dimension(:),allocatable,save::idea,idg,idcor,idcontr,indiceY
!   integer,dimension(:),allocatable,save::idlink,ntr
!   integer,dimension(:,:),allocatable,save::nmes
!   integer,dimension(:),allocatable,save::nvalSPL,nvalORD
!   double precision,dimension(:),allocatable,save :: seqMC
!   double precision,dimension(:),allocatable,save :: epsY
!   double precision,dimension(:,:),allocatable,save::zitr
!   double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2
!   double precision,dimension(:),allocatable,save::Tmm,Tmm1,&
!        Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02,Tmm03,Tim0, &
!        Tim01,Tim02,Tim03,Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1,&
!        Timt2,Timt3
!   double precision,dimension(:,:),allocatable,save::Tmm_st2,Tmm1_st2,&
!        Tmm2_st2,Tmm3_st2,Tmm0_st2,Tmm01_st2,Tmm02_st2,Tmm03_st2
!   integer,dimension(:),allocatable,save::fix
!   double precision,dimension(:),allocatable,save::bfix
!   integer,save::idst        
!   integer,dimension(2),save::nXcl        
!   integer,dimension(2),save::id_nXcl 
!   double precision,dimension(:,:),allocatable,save::Xcl_Ti,Xcl_GK,Xcl0_GK 
  
! end module mla1


subroutine loglik1(Y0,X0,Tentr0,Tevt0,Devt0,ind_survint0 &
     ,idea0,idg0,idcor0,idcontr0,idsurv0,idtdv0 &
     ,typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0 &
     ,ny0,ns0,nv0,nobs0,nea0,nmes0,idiag0,ncor0,nalea0&
     ,npm0,b0,nfix0,bfix0,epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0 &
     ,nvalSPLORD0,fix0,methInteg0,nMC0,dimMC0,seqMC0 &
     ,idst0,nXcl0,Xcl_Ti0,Xcl_GK0,loglik)

  use modirtsre

  IMPLICIT NONE

  !Declaration des variables en entree
  integer,intent(in)::nv0,ny0,nMC0,methInteg0,dimMC0,nfix0
  integer, intent(in)::ns0,nobs0,idiag0,npm0,nea0,ncor0,nalea0
  integer,intent(in)::idtrunc0,logspecif0,nbevt0
  double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
  integer, dimension(ns0),intent(in)::ind_survint0,Devt0
  integer, dimension(nv0),intent(in)::idtdv0,idsurv0
  integer,dimension(nbevt0),intent(in)::typrisq0,nz0
  double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0    
  double precision,dimension(ny0),intent(in)::epsY0
  integer, dimension(ny0),intent(in)::idlink0,nbzitr0,nvalSPLORD0
  double precision,dimension(maxval(nbzitr0),ny0),intent(in)::zitr0
  integer,dimension(nobs0),intent(in)::indiceY0
  double precision,dimension(sum(nvalSPLORD0(:))),intent(in)::uniqueY0
  integer, dimension(nv0),intent(in)::idea0,idg0,idcor0,idcontr0
  integer,dimension(ns0,ny0)::nmes0   
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  integer,dimension(npm0+nfix0),intent(in)::fix0
  double precision,dimension(dimMC0*nMC0),intent(in)::seqMC0
  integer,intent(in)::idst0   
  integer,dimension(2),intent(in)::nXcl0   
  double precision,dimension(ns0,nXcl0(1)),intent(in)::Xcl_Ti0 
  double precision,dimension(15*ns0,nXcl0(2)),intent(in)::Xcl_GK0
  double precision, dimension(npm0), intent(in) :: b0
  double precision, dimension(nfix0), intent(in) :: bfix0
  
  !Declaration des variables en sortie
  double precision,intent(out)::loglik

  !Variables locales
  integer::jtemp,i,j,npm,ier,k,ktemp,yk,k1,k2,mi,nbfix,p
  integer::ke,sumnrisq,it,npmtot0
  double precision::eps,ca,cb,dd
  double precision,external::vrais1


  maxmes=0
  do i=1,ns0
     mi=sum(nmes0(i,:))
     if (mi.gt.maxmes) then
        maxmes=mi
     end if
  end do


  allocate(rangeY(ny0),minY(ny0),maxY(ny0),idlink(ny0),ntr(ny0),epsY(ny0))

  nySPL=0
  nyORD=0
  rangeY=0
  epsY=epsY0
  do k=1,ny0
     idlink(k)=idlink0(k)
     minY(k)=zitr0(1,k)
     maxY(k)=zitr0(nbzitr0(k),k)
     if (idlink(k).eq.2) then
        nySPL=nySPL+1
     end if
     if (idlink(k).eq.3) then
        nyORD=nyORD+1
     end if
  end do
  
  if(nySPL.gt.0) then 
     allocate(nvalSPL(nySPL))
     nvalSPL=0
  else
     allocate(nvalSPL(1))
     nvalSPL(1) = 0
  end if

  if(nyORD.gt.0) then
     allocate(nvalORD(nyORD))
     nvalORD=0
  else
     allocate(nvalORD(1))
     nvalORD(1) = 0
  end if

  k1=0
  k2=0
  do k=1,ny0
     if(idlink(k).eq.2) then
        k1=k1+1
        nvalSPL(k1)=nvalSPLORD0(k)
     else if (idlink(k).eq.3) then
        k2=k2+1
        nvalORD(k2)=nvalSPLORD0(k)
     end if
  end do
  ntotvalSPL=sum(nvalSPL(:))
  ntotvalORD=sum(nvalORD(:))

  methInteg = methInteg0
  nMC = nMC0
  
  if(all(idlink.ne.2)) then
     allocate(zitr(1,1))
     allocate(mm(1),mm1(1),mm2(1),im(1),im1(1),im2(1))
     mm(1)=0.d0
     mm1(1)=0.d0
     mm2(1)=0.d0
     im(1)=0.d0
     im1(1)=0.d0
     im2(1)=0.d0
  else
     allocate(zitr(-1:(maxval(nbzitr0)+2),nySPL))
     allocate(mm(ntotvalSPL),mm1(ntotvalSPL),mm2(ntotvalSPL),im(ntotvalSPL),im1(ntotvalSPL),im2(ntotvalSPL))
  end if


  zitr=0.d0  
  k1=0
  k2=0
  do k=1,ny0
     if (idlink(k).eq.0) ntr(k)=2     
     if (idlink(k).eq.1) ntr(k)=4
     if (idlink(k).eq.2) then
        k1=k1+1
        ntr(k)=nbzitr0(k)+2

        zitr(1:nbzitr0(k),k1)=zitr0(1:nbzitr0(k),k)
        zitr(-1,k1)=zitr(1,k1)
        zitr(0,k1)=zitr(1,k1)
        zitr(ntr(k)-1,k1)=zitr(ntr(k)-2,k1)
        zitr(ntr(k),k1)=zitr(ntr(k)-1,k1)
     end if
     if (idlink(k).eq.3) then
        k2 = k2+1
        ntr(k) = nvalORD(k2)-1
     end if
  end do

  ntrtot = sum(ntr)
  allocate(Y(nobs0),X(nobs0,nv0),uniqueY(ntotvalSPL+ntotvalORD) &
       ,idea(nv0),idg(nv0),idcor(nv0),idcontr(nv0),nmes(ns0,ny0) &
       ,indiceY(nobs0))

  allocate(Tsurv0(ns0),Tsurv(ns0),Tsurvint(ns0),ind_survint(ns0),Devt(ns0))
  allocate(Tsurv0_st2(ns0,15),Tsurv_st2(ns0,15))
  allocate(typrisq(nbevt0),nz(nbevt0),nprisq(nbevt0),nevtparx(nv0),nxcurr(nv0))
  allocate(idsurv(nv0),idtdv(nv0))

  ! zi : contient noeuds pour hazard (ou min et max si Weibull)
  if(any(typrisq0.eq.3)) then
     allocate(zi(-2:maxval(nz0)+3,nbevt0))
  else
     allocate(zi(maxval(nz0),nbevt0))
  end if
  
  allocate(Xcl_Ti(ns0,nXcl0(1)),Xcl_GK(15*ns0,nXcl0(1)),Xcl0_GK(15*ns0,nXcl0(1)))

  eps=1.d-20

  ! enregistrement pour les modules
  nbevt=nbevt0    
  typrisq=typrisq0
  idtrunc=idtrunc0
  Tsurv0=Tentr0   
  Tsurv=Tevt0    
  devt=devt0    
  ind_survint=ind_survint0
  logspecif=logspecif0 
  Tsurvint=Tsurv
  ny=ny0
  ns=ns0
  nv=nv0
  nobs=nobs0
  ncor=ncor0
  nalea=nalea0
  idiag=idiag0
  idst=idst0  
  nXcl=nXcl0
  if(idst.eq.2) then
     Xcl_Ti=Xcl_Ti0
     Xcl_GK=Xcl_GK0(:,1:nXcl(1))
     if (idtrunc.eq.1) then
        Xcl0_GK=Xcl_GK0(:,(nXcl(1)+1):nXcl(2))
     end if
     do i=1,ns0  
        do p=1,15
           Tsurv_st2(i,p) = Xcl_GK(15*(i-1)+p,1)
           if (idtrunc.eq.1) then
              Tsurv0_st2(i,p) = Xcl0_GK(15*(i-1)+p,1) 
           end if
        end do
     end do
  end if
  
  
  !     if (verbose==1) write(*,*)'ntotvalSPL',ntotvalSPL

  if (ntotvalSPL+ntotvalORD.gt.0) uniqueY(1:ntotvalSPL+ntotvalORD)=uniqueY0(1:ntotvalSPL+ntotvalORD)

  nmes=0
  Y=0.d0
  X=0.d0
  idea=0
  idg=0
  idcor=0
  idcontr=0
  idsurv=0
  idtdv=0
  ktemp=0

  do k=1,nv
     idsurv(k)=idsurv0(k)
     idtdv(k)=idtdv0(k)
     idea(k)=idea0(k)
     idg(k)=idg0(k)
     idcor(k)=idcor0(k)
     idcontr(k)=idcontr0(k)

     jtemp=0
     DO i=1,ns
        do yk=1,ny            
           if (k.eq.1) then
              nmes(i,yk)=nmes0(i,yk)   !dim(nmes)=ns*ny    
              do j=1,nmes(i,yk)
                 jtemp=jtemp+1
                 Y(jtemp)=Y0(jtemp)
                 indiceY(jtemp)=indiceY0(jtemp)
                 ktemp=ktemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           else
              do j=1,nmes(i,yk)
                 ktemp=ktemp+1
                 jtemp=jtemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           end if
        end do
     end do
  end do
          
  id_nXcl(1)=sum(idg)-1 !sans intercept
  id_nXcl(2)=sum(idea)
  
  ! definition Tsurvint 
  nvdepsurv=0
  if(sum(ind_survint).gt.0) then
     nvdepsurv=1
     do k=1,nv
        if (idtdv(k).eq.1) then
           it=0
           do i=1,ns
              Tsurvint(i)=X(it+1,k)
              it=it+maxval(nmes(i,:))
           end do
        end if
     end do
  end if


  ! prm fixes
  npmtot0 = npm0+nfix0
  allocate(fix(npmtot0))
  fix=0
  fix(1:npmtot0)=fix0(1:npmtot0)
  nbfix=sum(fix)
  if(nbfix.eq.0) then
     allocate(bfix(1))
  else
     allocate(bfix(nbfix))
  end if
  bfix(1:nbfix)=bfix0(1:nbfix)

  ! creation des parametres

  nprisq=0
  nrisqtot=0

  do ke=1,nbevt

     nz(ke)=nz0(ke) ! nb de noeuds pour hazard (ou 2 pr Weibull)

     if (typrisq(ke).eq.1) then
        nprisq(ke)=nz(ke)-1
     end if
     if (typrisq(ke).eq.2) then
        nprisq(ke)=2
     end if
     if (typrisq(ke).eq.3) then
        nprisq(ke)=nz(ke)+2
     end if

     nrisqtot = nrisqtot+nprisq(ke)  ! nb total de prm pour hazards
     zi(1:nz(ke),ke)=zi0(1:nz(ke),ke)
  end do

  ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
  nxevt=0
  nevtparx=0
  do j=1,nv

     if(idtdv(j).ne.1) then

        if(idsurv(j).eq.1) then
           nevtparx(j) = 1
           nxevt = nxevt + 1
        end if
        if(idsurv(j).eq.2) then 
           nevtparx(j) = nbevt
           nxevt = nxevt + 1
        end if
     end if

  end do

  nvarxevt = sum(nevtparx) + nvdepsurv

  nea=0
  nef=0
  ncontr=0
  do k=1,nv
     if (idg(k).eq.1) then
        nef = nef + 1
     end if
     nea=nea+idea(k)
     ncontr=ncontr+idcontr(k)*(ny-1) 
  end do
  nef = nef - 1 !intercept pas estime

  nasso = nbevt*nea
  if (idst.eq.2) then
    nasso = nbevt
  end if

  if (idiag.eq.1) then
     nvc=nea-1
  else if(idiag.eq.0) then
     nvc=(nea+1)*nea/2-1
  end if

  npmtot = nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+ny

  ! points qmc
  if(methInteg.ne.3) then 
     allocate(seqMC(1))
  else
     allocate(seqMC(dimMC0*nMC))
     seqMC = seqMC0(1:dimMC0*nMC) 
  end if



  ! creer base de splines si au moins un hazard splines
  if(any(typrisq.eq.3)) then
     allocate(Tmm(ns*nbevt),Tmm1(ns*nbevt),Tmm2(ns*nbevt),Tmm3(ns*nbevt),        &
              Tim(ns*nbevt),Tim1(ns*nbevt),Tim2(ns*nbevt),Tim3(ns*nbevt),        &
              Tmm0(ns*nbevt),Tmm01(ns*nbevt),Tmm02(ns*nbevt),Tmm03(ns*nbevt),    &
              Tim0(ns*nbevt),Tim01(ns*nbevt),Tim02(ns*nbevt),Tim03(ns*nbevt),    &
              Tmmt(ns*nbevt),Tmmt1(ns*nbevt),Tmmt2(ns*nbevt),Tmmt3(ns*nbevt),    &
              Timt(ns*nbevt),Timt1(ns*nbevt),Timt2(ns*nbevt),Timt3(ns*nbevt),    &
              Tmm_st2(15,ns*nbevt),Tmm1_st2(15,ns*nbevt),Tmm2_st2(15,ns*nbevt),Tmm3_st2(15,ns*nbevt),        &
              Tmm0_st2(15,ns*nbevt),Tmm01_st2(15,ns*nbevt),Tmm02_st2(15,ns*nbevt),Tmm03_st2(15,ns*nbevt))

     do ke=1,nbevt
        if(typrisq(ke).eq.3) then
           call splines_irtsre(ke)
           call splines_irtsre2(ke)
        end if
     end do

  end if



  ! base de splines transfos
  if (any(idlink.eq.2)) then 
     call design_splines_irtsre(ier)
     if (ier.eq.-1) then
        loglik=-1.d9
        go to 1589
     end if
  end if

  ! calcul de la vraisemblance
  loglik = vrais1(b0,npm0)

  
1589 continue

  if (any(typrisq.eq.3)) then
     deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
          Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
          Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3,    &
          Tmm_st2,Tmm1_st2,Tmm2_st2,Tmm3_st2,    &   
          Tmm0_st2,Tmm01_st2,Tmm02_st2,Tmm03_st2)
  endif

  deallocate(Tsurv0,Tsurv,Tsurvint &
       ,Tsurv0_st2, Tsurv_st2 & 
       ,ind_survint,zi,devt,typrisq,nz,nprisq,idsurv,idtdv &
       ,nevtparx,nxcurr)

  deallocate(Y,X,idea,idg,idcor,idcontr,nmes,uniqueY,indiceY,ntr)


  deallocate(zitr,mm,mm1,mm2,im,im1,im2,minY,maxY,rangeY,idlink,nvalSPL,nvalORD,epsY)

  deallocate(fix,bfix,seqMC)

  deallocate(Xcl_Ti,Xcl_GK,Xcl0_GK) 
  
  return
  
end subroutine loglik1







!-----------------------------------------------------------
!                       VRAIS1_i
!------------------------------------------------------------


double precision function vrais1_i(b,npm,i) 

  use modirtsre
  use optim

  IMPLICIT NONE
  integer ::i,j,k,l,m,jj,npm,ll,ii,numSPL,ykord
  integer ::ier,kk,j1,j2,sumMesYk,yk,sumntr,ke,sumnrisq
  integer::nevtxcurr,nxevtcurr
  double precision,dimension(maxmes,nv) ::X00
  double precision,dimension(maxmes,nea) ::Z
  double precision,dimension(maxmes,(ncontr+sum(idcontr)))::X01
  double precision,dimension(ncontr+sum(idcontr))::b01
  double precision,dimension(nea,nea) ::Ut
  double precision,dimension(maxmes,maxmes) ::VC,Corr
  double precision,dimension(npm) :: b
  double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
  double precision,dimension(nv) :: b0
  double precision,dimension(npmtot)::b1
  double precision,dimension(nxevt)::Xevt,bevt
  double precision,dimension(nbevt)::bevtint
  double precision,dimension(maxval(nprisq))::brisq
  double precision::basso 
  double precision,dimension(nef)::beta_ef 

  double precision :: eps,det,som,eta0
  double precision ::Y4,jacobien,beta_densite,ytemp
  double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
  double precision,dimension(-1:maxval(ntr)-3)::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,cc1
  double precision,dimension(nea)::ui,usim
  double precision,dimension(maxmes)::wi,wsim
  double precision,dimension(ny)::asim
  double precision::ai,binf,bsup
  double precision,dimension(nbevt)::risq,surv,surv0,survint
  double precision::SX,x22,div,vrais_Y,vrais_surv,varexpsurv
  double precision::surv0_glob,surv_glob,fevt,easurv
  double precision,external::alnorm
  double precision::pred_cl_Ti 
  double precision::som_T0,som_Ti
  
  ! definir le nombre total de mesures pour le sujet i : nmestot (valable que pour cette fonction)

  ! if (verbose==1) write(*,*)'i',i 
  b1=0.d0
  eps=1.D-20
  l=0
  m=0
  do k=1,npmtot
     if(fix(k).eq.0) then
        l=l+1
        b1(k)=b(l)
     end if
     if(fix(k).eq.1) then
        m=m+1
        b1(k)=bfix(m)
     end if
  end do

  !----------- rappel des parametres utilises ---------



  !        write(*,*)'i',i,nmescur,nmescur + nmes(i)


  Ut=0.d0
  Ut(1,1)=1.d0
  if (nea>1) then 

     If (idiag.eq.1) then
        do j=2,nea
           do k=2,nea
              if (j.eq.k) then
                 Ut(j,k)=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+j-1)
              else
                 Ut(j,k)=0.d0
              end if
           end do
        end do
     end if

     If (idiag.eq.0) then
        do j=2,nea
           do k=1,j
              Ut(j,k)=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+k-1+j*(j-1)/2)
           end do
        end do
     end if

  end if
  
  vrais_Y=0.d0
  jacobien=0.d0
  ! -------- creation de Vi = ZiGZi'+se*seIni ----------
  ! creation de Zi

  Z=0.d0
  l=0
  do k=1,nv
     if (idea(k).eq.1) then
        l=l+1
        do j=1,sum(nmes(i,:))
           Z(j,l)=dble(X(nmescur+j,k))  
        end do
     end if
  end do

  !matrice Corr variance de BM/AR

  Corr=0.d0
  tcor=0.d0
  if (ncor.gt.0) then
     do k=1,nv
        if (idcor(k).eq.1) then
           do j=1,sum(nmes(i,:))
              tcor(j) = X(nmescur+j,k)
           end do
        end if
     end do
     do j1=1,sum(nmes(i,:))
        do j2=1,sum(nmes(i,:))
           if (ncor.eq.1) then 
              Corr(j1,j2) = Corr(j1,j2)+b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                   b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)*min(tcor(j1),tcor(j2))
           else if (ncor.eq.2) then
              Corr(j1,j2) = Corr(j1,j2)+b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                   b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                   exp(-b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1)*abs(tcor(j1)-tcor(j2)))
           end if
        end do
     end do

     ! passer en cholesky si on a de l ordinal
     if(any(idlink.eq.3)) then
        jj=0
        Vi=0.d0
        do j=1,sum(nmes(i,:))
           do k=j,sum(nmes(i,:))
              jj=j+k*(k-1)/2
              Vi(jj)=Corr(j,k)
           end do
        end do

        CALL DMFSD(Vi,sum(nmes(i,:)),EPS,IER)

        Corr=0.d0
        do j=1,sum(nmes(i,:))
           do k=1,j
              Corr(j,k)=Vi(k+j*(j-1)/2)
           end do
        end do
     end if
  end if

  ! creation de Y1
  Y1=0.d0
  splaa=0.d0

  sumMesYk = 0
  sumntr=0
  numSPL=0
  do yk=1,ny

     if (idlink(yk).eq.0) then  ! Linear link

        do j=1,nmes(i,yk)
           Y1(sumMesYk+j)=(dble(Y(nmescur+sumMesYk+j))-b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)) &
                /abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))
            
           jacobien = jacobien - log(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))
        end do

     else if (idlink(yk).eq.1) then  ! Beta link


        aa1=exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1))/ &
             (1+exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)))
        bb1=exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))/ &
             (1+exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2)))
        bb1=aa1*(1.d0-aa1)*bb1

        cc1=abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+3))

        dd1=abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+4))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1

        do j=1,nmes(i,yk)

           ytemp=(dble(Y(nmescur+sumMesYk+j))-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
           Y1(sumMesYk+j)=(betai(aa,bb,ytemp)-cc1)/dd1


           if (Y1(sumMesYk+j).eq.999.d0) then
              vrais1_i=-1.d9
              !print*,"-1.d9 Y1=999"
              goto 654
           end if

           jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
           jacobien=jacobien-log(abs(maxY(yk)-minY(yk)+2*epsY(yk)))
        end do

     else if (idlink(yk).eq.2) then ! Splines link
        numSPL=numSPL+1

        splaa=0.d0
        eta0=0.d0
        eta0=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)

        do kk=2,ntr(yk)
           splaa(kk-3)=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+kk)&
                *b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+kk)
        end do
        !if(i==1 .and. id==0 .and. jd==0) print*,"eta0=",eta0,"splaa=",sqrt(splaa)
        do j=1,nmes(i,yk)
           ll=0
           !if(i==1 .and. id==0 .and. jd==0) print*,"Y=",Y(nmescur+sumMesYk+j)
           if (Y(nmescur+sumMesYk+j).eq.zitr(ntr(yk)-2,numSPL)) then
              ll=ntr(yk)-3
           end if

           som=0.d0
           do kk = 2,ntr(yk)-2
              if ((Y(nmescur+sumMesYk+j).ge.zitr(kk-1,numSPL)).and. &
                   (Y(nmescur+sumMesYk+j).lt.zitr(kk,numSPL))) then
                 ll=kk-1
              end if
           end do

           if (ll.lt.1.or.ll.gt.ntr(yk)-3) then          
              vrais1_i=-1.d9
              print*,"-1.d9 ll<1 ou ll>ntr-3",ll!," ntr=",ntr(yk)," numSPL=",numSPL," y=",Y(nmescur+sumMesYk+j)
              goto 654
           end if
           if (ll.gt.1) then
              do ii=2,ll
                 som=som+splaa(ii-3)
              end do
           end if



           Y1(sumMesYk+j)=eta0+som +splaa(ll-2)*im2(indiceY(nmescur+sumMesYk+j)) &
                +splaa(ll-1)*im1(indiceY(nmescur+sumMesYk+j))&
                + splaa(ll)*im(indiceY(nmescur+sumMesYk+j))

           jacobien = jacobien + log(splaa(ll-2)*mm2(indiceY(nmescur+sumMesYk+j)) &
                +splaa(ll-1)*mm1(indiceY(nmescur+sumMesYk+j))&
                +splaa(ll)*mm(indiceY(nmescur+sumMesYk+j)))

           !print*,"jac=",jacobien
           !print*,"ll =",ll
           !print*,"splaa =",splaa(ll-2)
           !print*,"nmescur+sumMesYk+j =",nmescur+sumMesYk+j
           !print*,"indiceY= ",indiceY(nmescur+sumMesYk+j)
           !print*,"mm =",mm(nmescur+sumMesYk+j)
           !print*,"mm1 =",mm1(nmescur+sumMesYk+j)
           !print*,"mm2 =",mm2(nmescur+sumMesYk+j)
           !print*,"Y =",Y(nmescur+sumMesYk+j)

           !                write(*,*)'Y',Y1(sumMesYk+j),sumMesYk,yk,j,jacobien
        end do
     else if (idlink(yk).eq.3) then  ! Threshold link
        do j=1,nmes(i,yk)
           Y1(sumMesYk+j)=Y(nmescur+sumMesYk+j)
           !if(nmescur.lt.15) then
           !   print*,"nmescur=",nmescur," sumMesYk=",sumMesYk," j=",j
           !   print*,"Y=",Y(nmescur+sumMesYk+j), "  Y1=",Y1(sumMesYk+j)
           !end if
        end do
     end if
     sumMesYk=sumMesYk+nmes(i,yk)
     sumntr=sumntr+ntr(yk)
     
  end do !fin boucle yk


  !         if (i.lt.3)then
  !            write(*,*)'nmes',nmes(i),b1((nef+ncontr+nvc+1):npm),nef+ncontr
  !            write(*,*)'Y1',Y1
  !         end if




  ! contribution individuelle a la vraisemblance
  ! print*,"i=",i," -ni*log(2pi)=",-sum(nmes(i,:))*dlog(dble(2*3.14159265)), " log(det)=",det
  ! print*,"Vi=",VC
  ! sans classes latentes : ng=1
  vrais1_i=0.d0


  b0=0.d0
  b01=0.d0
  l=0
  m=0
  X00=0.d0
  X01=0.d0
  do k=1,nv
     if (idg(k).ne.0) then
        l=l+1
        do j=1,sum(nmes(i,:))
           X00(j,l)=dble(X(nmescur+j,k))  
        end do
        ! idg ne 0 pour l'intercept forcement donc on met le parm a 0
        if (k.eq.1) then
           b0(l)=0.d0
        else
           b0(l)=b1(nrisqtot+nvarxevt+nasso+l-1)
        end if
     end if

     !contrast : 
     if (idcontr(k).ne.0) then
        m=m+1
        sumMesYk=0
        do yk=1,ny
           ! creation matrice design des contrastes: X01
           do j=1,nmes(i,yk)
              X01(sumMesYk+j,(m-1)*ny+yk) = dble(X(nmescur+sumMesYk+j,k))
           end do
           sumMesYk=sumMesYk+nmes(i,yk)
           ! creation vecteur parms des contrastes: b01
           if (yk<ny) THEN
              b01((m-1)*ny+yk)=b1(nrisqtot+nvarxevt+nasso+nef+(m-1)*(ny-1)+yk)
           else
              b01((m-1)*ny+ny) =-sum(b1(nrisqtot+nvarxevt+nasso+nef+(m-1)*(ny-1)+1 &
                   :nef+(m-1)*(ny-1)+ny-1))
           end if
        end do
     end if
  end do

  som=0.d0
  som_T0=0.d0 !delayed entry
  som_Ti=0.d0
  do l=1,nMC

     vrais_Y=1.d0
     mu=0.d0

!!!!!!!!! MC pour EA et BM/AR !!!!!!!!!

     if(methInteg.eq.1) then 
        ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!

        ! simuler les effets aleatoires
        if(nea.gt.0) then
           x22=0.d0
           SX=1.d0
           do j=1,nea
              call bgos(SX,0,usim(j),x22,0.d0)
              !print*,"usim=",usim(j)
           end do
           ui=0.d0
           ui=matmul(Ut,usim)
           !print*,"usim=",usim(j)," ui=",ui, " Ut=",Ut, "  nea=",nea
        end if

        ! simuler le BM ou AR
        if(ncor.gt.0) then
           x22=0.d0
           SX=1.d0
           do j=1,sum(nmes(i,:))
              call bgos(SX,0,wsim(j),x22,0.d0)
           end do
           wi=0.d0
           wi=matmul(Corr,wsim)
        end if

     else if(methInteg.eq.2) then 
        ! !!!!!!!!!!!!! MCA !!!!!!!!!!!!!


        if(mod(l,2).eq.0) then
           ! si l est pair on prend l'oppose des precedents
           ui = -ui
           wi = -wi
        else
           ! sinon on simule des nouveaux

           ! simuler les effets aleatoires
           if(nea.gt.0) then
              x22=0.d0
              SX=1.d0
              do j=1,nea
                 call bgos(SX,0,usim(j),x22,0.d0)
              end do
              ui=0.d0
              ui=matmul(Ut,usim)
           end if

           ! simuler le BM ou AR
           if(ncor.gt.0) then
              x22=0.d0
              SX=1.d0
              do j=1,sum(nmes(i,:))
                 call bgos(SX,0,wsim(j),x22,0.d0)
              end do
              wi=0.d0
              wi=matmul(Corr,wsim)
           end if

        end if

     else 
        ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!

        ! simuler les effets aleatoires
        if(nea.gt.0) then
           usim=0.d0
           do j=1,nea
              usim(j)=seqMC(nMC*(j-1)+l)
           end do
           ui=0.d0
           ui=matmul(Ut,usim)
        end if

        ! simuler le BM ou AR
        if(ncor.gt.0) then
           wsim=0.d0
           do j=1,sum(nmes(i,:))
              wsim(j)=seqMC(nMC*(nea+j-1)+l)
           end do
           wi=0.d0
           wi=matmul(Corr,wsim)
        end if


     end if ! fin if methInteg


     ! esperance conditionnelle
     mu = matmul(X00,b0)+matmul(X01,b01)+matmul(Z,ui)

     if(ncor.gt.0) mu = mu+wi
     !if(i.lt.4) then
     !print*,"i=",i," nmes=",nmes(i,1)," nmescur=",nmescur
     !print*,"Xb=",matmul(X00,b0)+matmul(X01,b01)
     !print*,"mu=",mu
     !print*,"ui=",ui, " wi=",wi," xb=",matmul(X00,b0)+matmul(X01,b01)
     !end if
     sumMesYk=0
     sumntr=0
     ykord=0
     do yk =1,ny

        if(idlink(yk).eq.3) then
           !! yk est ordinal
           ykord = ykord + 1

           !! MC pour simuler l'EA specifique au test
           ai=0.d0
           if(nalea.gt.0) then
              if(methInteg.eq.1) then
                 !! MCO
                 call bgos(SX,0,asim(yk),x22,0.d0)
                 ai = abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
              else if(methInteg.eq.2) then
                 !! MCA
                 if(mod(l,2).eq.0) then
                    ! si l est pair on prend l'oppose du precedent
                    ai = -abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
                 else
                    call bgos(SX,0,asim(yk),x22,0.d0)
                    ai = abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
                 end if
              else
                 !! QMC
                 asim(yk) = seqMC(nMC*(nea+sum(nmes(i,:)))+l)
                 ai = abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
              end if
           end if
!if(i.lt.4) print*,"i=",i," avant do j, vrais_Y=",vrais_Y
           do j=1,nmes(i,yk)

              !! on ajoute ai a mu
              if(nalea.gt.0) mu(sumMesYk+j) = mu(sumMesYk+j)+ai

              !! trouver binf et bsup tq binf < lambda + epsilon < bsup

              !! initialiser au premier seuil
              binf = b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)
              bsup = binf

              !! si Y>minY ajouter b1(..)^2
              if(indiceY(nmescur+sumMesYk+j).gt.1) then
                 do ll=2,min(indiceY(nmescur+sumMesYk+j),ntr(yk))
                    bsup = bsup + b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+ll)**2
                    if(ll.lt.indiceY(nmescur+sumMesYk+j)) then
                       binf = binf + b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+ll)**2
                    end if
                 end do
              end if
              !if(i.lt.4)print*,"y=",Y1(sumMesYk+j)," indiceY=",indiceY(nmescur+sumMesYk+j), " Y=",Y(nmescur+sumMesYk+j)," mu=",mu(sumMesYk+j)
              !print*," binf=",binf," bsup=",bsup
              !! centrer et standardiser
              binf = (binf - mu(sumMesYk+j))/abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk))
              bsup = (bsup - mu(sumMesYk+j))/abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk))

              !if(i.lt.4) print*,"j=",j," nvalORD=",nvalORD(ykord)," binf=",binf," bsup=",bsup
              
              if(indiceY(nmescur+sumMesYk+j).eq.1) then
                 !! si Y=minY
                 vrais_Y = vrais_Y * alnorm(binf,.false.)
              else if(indiceY(nmescur+sumMesYk+j).eq.nvalORD(ykord)) then
                 !! si Y=maxY
                 vrais_Y = vrais_Y * (1.d0-alnorm(bsup,.false.))
              else
                 !! minY < Y < maxY
                 vrais_Y = vrais_Y * (alnorm(bsup,.false.)-alnorm(binf,.false.))
              end if
              !if(i.lt.4) print*,"vrais_Y=",vrais_Y

           end do


        else
           
           !! yk est continu
           !print*,"MC continu"
           if(nmes(i,yk).gt.0) then
              !! variance de Y|ui,wi
              VC=0.d0
              do j1=1,nmes(i,yk)
                 VC(j1,j1) = b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk)**2 !variance de l'erreur yk
                 if (nalea.eq.ny) then ! intercept aleatoire de yk
                    do j2=1,nmes(i,yk)
                       VC(j1,j2) = VC(j1,j2) + &
                            b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk)**2
                    end do
                 end if
              end do
              
              ! Vi en vecteur
              jj=0
              Vi=0.d0
              do j1=1,nmes(i,yk)
                 do j2=j1,nmes(i,yk)
                    jj=j1+j2*(j2-1)/2
                    Vi(jj)=VC(j1,j2)
                 end do
              end do
              
              ! inversion
              CALL dsinv(Vi,nmes(i,yk),eps,ier,det)
              if (ier.eq.-1) then
                 vrais1_i=-1.d9
                 print*,"-1.d9 dsinv continu MC"
                 !print*,"b=",b
                 !print*,"bfix=",bfix
                 !print*,"fix=",fix
                 goto 654
              end if
              
              ! retransformation du vecteur Vi en matrice :
              VC=0.d0
              do j1=1,nmes(i,yk)
                 do j2=1,nmes(i,yk)
                    if (j2.ge.j1) then
                       VC(j1,j2)=Vi(j1+j2*(j2-1)/2)
                    else
                       VC(j1,j2)=Vi(j2+j1*(j1-1)/2)
                    end if
                 end do
              end do
           
              ! calcul de la vrais
              Y2=0.d0
              Y3=0.d0
              Y4=0.d0
              do j=1,nmes(i,yk)
                 Y2(j) = Y1(sumMesYk+j)-mu(sumMesYk+j)
              end do
              Y3=matmul(VC,Y2)
              Y4=DOT_PRODUCT(Y2,Y3)
              
              div = (dble(2*3.14159265)**(dble(nmes(i,yk))/2))*sqrt(exp(det))
              vrais_Y = vrais_Y * exp(-Y4/2.d0)/div

           end if
        end if
        
        sumMesYk = sumMesYk+nmes(i,yk)
        sumntr = sumntr + ntr(yk)
     end do ! fin boucle yk

!if(i.lt.4) print*,"avant survie, vrais_Y",vrais_Y
     ! partie survie

     if (nbevt.ne.0) then

        ! calcul de brisq en chaque composante et risq, surv et surv0 pour chaque evt
        risq=0.d0
        surv=0.d0
        surv0=0.d0
        survint=0.d0

        sumnrisq=0
        do ke=1,nbevt

           brisq=0.d0
           if (logspecif.eq.1) then
              do k=1,nprisq(ke)
                 brisq(k)=exp(b1(sumnrisq+k))
              end do
           else
              do k=1,nprisq(ke)
                 brisq(k)=b1(sumnrisq+k)*b1(sumnrisq+k)
              end do
           end if
           
           basso=0.d0 
           basso=b1(nrisqtot+nvarxevt+ke) !basso = eta = prm estime des EAs partages
           
           beta_ef=0.d0 
           beta_ef=b1(nrisqtot+nvarxevt+ke) !beta_ef = prms des EFs necessaires pr pred curlev
           do k=1,nef
              beta_ef(k) = b1(nrisqtot+nvarxevt+nasso+k)
           end do

           
           if (idst.eq.1) then  !bi
              call fct_risq_irtsre(i,ke,brisq,risq,surv,surv0,survint)  
           else if(idst.eq.2) then   !niv.courant du processus latent
              call fct_risq_irtsre_2(i,ke,brisq,basso,beta_ef,ui,risq,surv,surv0)  
           end if

           sumnrisq = sumnrisq + nprisq(ke)
        end do

 ! print*,"fct_risq ok"
        ! variables explicatives de la survie
        Xevt=0.d0
        bevt=0.d0
        bevtint=0.d0
        if (nxevt.ne.0) then    ! si varexpsurv

           m=0
           do ke=1,nbevt
              nevtxcurr=0
              do k=1,nv 

                 if (idtdv(k).ne.1) then

                    if (idsurv(k).eq.1) then  
                       m=m+1
                       bevt(m)=b1(nrisqtot+nevtxcurr+1)
                       Xevt(m)=X(nmescur+1,k)
                    else
                       if (idsurv(k).eq.2) then   
                          m=m+1
                          bevt(m)=b1(nrisqtot+nevtxcurr+ke)
                          Xevt(m)=X(nmescur+1,k)
                       end if
                    end if

                 else ! i.e timedepvar

                    if (idsurv(k).eq.1) then  
                       bevtint(ke)=b1(nrisqtot+nevtxcurr+1)
                    else 
                       if (idsurv(k).eq.2) then
                          bevtint(ke)=b1(nrisqtot+nevtxcurr+ke)
                       end if
                    end if
                 end if
                 nevtxcurr=nevtxcurr+nevtparx(k)
              end do
           end do

        end if

        ! vrais survie
        surv_glob=0.d0
        surv0_glob=0.d0
        varexpsurv=0.d0
        nxevtcurr=0
        fevt=0.d0
        pred_cl_Ti=0.d0
        m=0
        do ke=1,nbevt
        
           ! calculer Xevt * bevt
           varexpsurv=0.d0
           if (nxevt.ne.0) then   !si varexpsurv
              varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevt))&
                   ,bevt((nxevtcurr+1):(nxevtcurr+nxevt)))
           end if
           
           ! effets aleatoires partages 
           if (idst.eq.1) then
            easurv=0.d0
            if(nea.gt.0) then      !si EA
              easurv=DOT_PRODUCT(ui,b1((nrisqtot+nvarxevt+m+1)&
                :(nrisqtot+nvarxevt+m+nea)))
              m = m+nea
            end if
           end if
           
           ! avoir evt au temps Ti si Devt=1
           if (Devt(i).eq.ke) then     !si sujet i a evt ke
              if (idst.eq.1) then
                fevt=risq(ke)*exp(varexpsurv+easurv)   !fct de risq         
              else if (idst.eq.2) then
                ! niv courant process latent au tps Ti
                do ll=1,id_nXcl(1)
                  pred_cl_Ti = pred_cl_Ti + Xcl_Ti(i,1+ll) * beta_EF(ll)  ! X(t) %*% beta
                end do
                do ll=1,id_nXcl(2)
                  pred_cl_Ti = pred_cl_Ti + Xcl_Ti(i,1+nef+ll) * ui(ll)  ! Z(t) %*% ui
                end do
                
                pred_cl_Ti = EXP(pred_cl_Ti*b1(nrisqtot+nvarxevt+ke) )  ! exp(predcl*basso)
                fevt=risq(ke)*exp(varexpsurv)*pred_cl_Ti
              end if
              if (ind_survint(i).eq.1) then
                 fevt=fevt*exp(bevtint(ke))
              end if
           end if
           
           ! risque cumule jusque Ti
           if (idst.eq.1) then
            Surv_glob=surv_glob + survint(ke)*exp(varexpsurv+easurv) + &
                  exp(bevtint(ke)+varexpsurv+easurv)*(surv(ke)-survint(ke))     
           else if (idst.eq.2) then
            Surv_glob=surv_glob + surv(ke)*exp(varexpsurv)
           end if
           
           ! troncature : risque cumule au temps T0
           if (idtrunc.eq.1) then
            if (idst.eq.1) then
              surv0_glob=surv0_glob+surv0(ke)*exp(varexpsurv+easurv)    
            else if (idst.eq.2) then
              surv0_glob=surv0_glob+surv0(ke)*exp(varexpsurv)  
            end if
           end if
           
           nxevtcurr=nxevtcurr+nxevt
        end do

        ! vraisemblance de la partie survie
        vrais_surv = exp(-Surv_glob)

 ! print*,"vrais_surv ok"        
        if(Devt(i).gt.0) vrais_surv = vrais_surv * fevt

        if (idtrunc.eq.1) then
           !vrais_surv = vrais_surv / exp(-surv0_glob)
           som_T0 = som_T0 + exp(-surv0_glob)   !delayed entry
        end if
        
        som_Ti = som_Ti + vrais_surv
        
        ! vrais totale 
        som = som + vrais_Y * vrais_surv

     else !  pas de survie

        som = som + vrais_Y

     end if
     
     !if(l.lt.4) print*,"l=", l, " som=",som
 !    print*,"fin l=",l
  end do ! fin boucle nMC
  
  vrais1_i = vrais1_i + log(som) - log(dble(nMC)) + jacobien
  
  if (idtrunc.eq.1) then    !delayedentry
      vrais1_i = vrais1_i - log(som_T0) + log(dble(nMC))
  end if  
  
!  print*,"trunc ok"
  !print*,"i=",i,"som_Ti=",som_Ti,"som_T0=",som_T0," vrais1_i=",vrais1_i
  
  !print*,"i=",i,"som=",som," vrais_Y=",vrais_Y,"jac=",jacobien," vrais_surv=",vrais_surv," vrais1_i=",vrais1_i
  654 continue

  return

end function vrais1_i





double precision function vrais1(b,m)


  use modirtsre,only:ns,nmes,nmescur

  implicit none

  integer::m,i
  double precision::vrais1_i,temp
  double precision,dimension(m)::b
  

  nmescur=0
  vrais1=0.d0
  do i=1,ns
  
!     print*,"## NEW SUBJECT i=",i 
  
     temp=vrais1_i(b,m,i)  
     
     vrais1 = vrais1 + temp
     if (temp.eq.-1.d9 .or. temp/temp.ne.1) then 
        !if (temp/temp.ne.1) write(*,*)"i=",i,"vrais= ",temp
        !if (temp.eq.-1.d9) then 
        vrais1 = -1.d9
        !print*,"dans vrais i=",i," vrais1=",vrais1," m=",m," b=",b
        ! if(verbose==1) write(*,*)"i=",i,"vrais= ",temp
        goto 541
     end if
     nmescur = nmescur + sum(nmes(i,:))
  end do
  
541 continue
  return

end function vrais1
