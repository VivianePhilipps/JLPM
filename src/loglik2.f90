

subroutine loglik2(Y,X,Tsurv0,Tsurv,Devt,Tsurvint, &
     idea,idg,idcor,idcontr,idsurv,idtdv, &
     typrisq,nz,zi,nbevt,idtrunc,logspecif, &
     ny,ns,nv,nobs,nea,nmes,idiag,&
     npm,b,nfix,bfix,epsY,idlink,nbzitr,zitr,indiceY, &
     fix,methInteg,nMC,dimMC,seqMC, &
     idst,nXcl,Xcl_Ti,Xcl_GK, & !43 param jusqu ici
     maxmes,nparam,nxevt,nprisq,ntr,&
     nyORD,nvalORD,ntotvalSPL,Tmm,Tim,Tim0,Timt,mm,im, &
     Tmm_st2,Tmm0_st2,loglik) !61

  use optim

  IMPLICIT NONE

  !Declaration des variables en entree
  integer,intent(in)::nv,ny,nMC,methInteg,dimMC,nfix,maxmes
  integer, intent(in)::ns,nobs,idiag,npm,nea,nyORD,ntotvalSPL
  integer,intent(in)::idtrunc,logspecif,nbevt,nxevt
  integer,dimension(nyORD)::nvalORD
  !integer,dimension(nySPL)::nvalSPL
  double precision, dimension(ns),intent(in)::Tsurv0,Tsurv,Tsurvint
  integer, dimension(ns),intent(in)::Devt
  integer, dimension(nv),intent(in)::idtdv,idsurv
  integer,dimension(nbevt),intent(in)::typrisq,nz,nprisq
  double precision,dimension(maxval(nz),nbevt),intent(in)::zi    
  double precision,dimension(ny),intent(in)::epsY
!  double precision,dimension(2*ny),intent(in)::rangeY
  integer, dimension(ny),intent(in)::idlink,nbzitr,ntr
  double precision,dimension(maxval(nbzitr),ny),intent(in)::zitr
  integer,dimension(nobs),intent(in)::indiceY
 ! double precision,dimension(sum(nvalSPLORD(:))),intent(in)::uniqueY
  integer, dimension(nv),intent(in)::idea,idg,idcor,idcontr
  integer,dimension(ns,ny)::nmes   
  double precision,dimension(nobs),intent(in)::Y
  double precision,dimension(nobs,nv),intent(in)::X
  integer,dimension(npm+nfix),intent(in)::fix
  double precision,dimension(dimMC*nMC),intent(in)::seqMC
  integer,intent(in)::idst   
  integer,dimension(2),intent(in)::nXcl   
  double precision,dimension(ns,nXcl(1)),intent(in)::Xcl_Ti 
  double precision,dimension(15*ns,nXcl(2)),intent(in)::Xcl_GK
  double precision, dimension(npm), intent(in) :: b
  double precision, dimension(nfix), intent(in) :: bfix
  integer, dimension(10)::nparam
  double precision, dimension(3*ntotvalSPL)::mm,im
  double precision, dimension(4*ns*nbevt)::Tmm,Tim,Tim0,Timt
  double precision,dimension(4*15,ns*nbevt)::Tmm_st2,Tmm0_st2


  !Declaration des variables en sortie
  double precision,intent(out)::loglik

  !Variables locales
  integer::i,j,ier,k,yk,nbfix
  integer::ke,sumnrisq,npmtot,nvdepsurv
  integer::nrisqtot,nvarxevt,nasso,nef,ncontr,nvc,ncor,nalea,ntrtot
  double precision::eps,vrais2_i
  integer::nmescur,l
  integer ::m,jj,ll,ii,ykord
  integer ::kk,j1,j2,sumMesYk,sumntr
  integer::nevtxcurr,nxevtcurr
  double precision,dimension(maxmes,nparam(4)+1) ::X00
  double precision,dimension(nparam(4)+1) :: b0
  double precision,dimension(maxmes,nea) ::Z
  double precision,dimension(nea,nea)::Ut
  double precision,dimension(maxmes,ny*sum(idcontr))::X01
  double precision,dimension(ny*sum(idcontr))::b01
  double precision,dimension(maxmes,maxmes) ::VC,Corr
  double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
  double precision,dimension(nxevt)::Xevt,bevt
  double precision,dimension(nbevt)::bevtint
  double precision,dimension(maxval(nprisq))::brisq
  double precision::basso 
  double precision,dimension(nparam(4))::beta_ef 

  double precision :: det,som,eta0
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
  double precision,dimension(npm+nfix)::btot


  nrisqtot = nparam(1)
  nvarxevt = nparam(2)
  nasso = nparam(3)
  nef = nparam(4)
  ncontr = nparam(5)
  nvc = nparam(6)
  ncor = nparam(7)
  ntrtot = nparam(8)
  nalea = nparam(9)
  !ny = nparam(10)

  nvdepsurv = sum(idtdv)


  npmtot = npm+nfix
  btot=0.d0
  jj=0
  nbfix=0
  do j=1,npmtot
     if(fix(j).eq.0) then
        jj = jj+1
        btot(j) = b(jj)
     else
        nbfix = nbfix+1
        btot(j) = bfix(nbfix)
     end if
  end do

  Ut=0.d0
  Ut(1,1)=1.d0
  if (nea>1) then 

     If (idiag.eq.1) then
        do j=2,nea
           do k=2,nea
              if (j.eq.k) then
                 Ut(j,k)=btot(nrisqtot+nvarxevt+nasso+nef+ncontr+j-1)
              else
                 Ut(j,k)=0.d0
              end if
           end do
        end do
     end if

     If (idiag.eq.0) then
        do j=2,nea
           do k=1,j
              Ut(j,k)=btot(nrisqtot+nvarxevt+nasso+nef+ncontr+k-1+j*(j-1)/2)
           end do
        end do
     end if

  end if

  vrais_Y=0.d0
  jacobien=0.d0

  nmescur=0
  loglik=0.d0
  do i=1,ns

          !print*,"## NEW SUBJECT i=",i 


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
                 Corr(j1,j2) = Corr(j1,j2)+btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                      btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)*min(tcor(j1),tcor(j2))
              else if (ncor.eq.2) then
                 Corr(j1,j2) = Corr(j1,j2)+btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                      btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                      exp(-btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1)*abs(tcor(j1)-tcor(j2)))
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
     jacobien=0.d0
     splaa=0.d0

     sumMesYk = 0
     sumntr=0
     do yk=1,ny

        if (idlink(yk).eq.0) then  ! Linear link

           do j=1,nmes(i,yk)
              Y1(sumMesYk+j)=(dble(Y(nmescur+sumMesYk+j))-btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)) &
                   /abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))

              jacobien = jacobien - log(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))
           end do

        else if (idlink(yk).eq.1) then  ! Beta link


           aa1=exp(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1))/ &
                (1+exp(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)))
           bb1=exp(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))/ &
                (1+exp(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2)))
           bb1=aa1*(1.d0-aa1)*bb1

           cc1=abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+3))

           dd1=abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+4))

           aa=aa1*aa1*(1-aa1)/bb1-aa1
           bb=aa*(1-aa1)/aa1

           do j=1,nmes(i,yk)

              ytemp=(dble(Y(nmescur+sumMesYk+j))-zitr(1,yk)+epsY(yk))/(zitr(2,yk)-zitr(1,yk)+2*epsY(yk))
              Y1(sumMesYk+j)=(betai(aa,bb,ytemp)-cc1)/dd1


              if (Y1(sumMesYk+j).eq.999.d0) then
                 loglik=-1.d9
                 !print*,"-1.d9 Y1=999"
                 goto 654
              end if

              jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
              jacobien=jacobien-log(abs(zitr(2,yk)-zitr(1,yk)+2*epsY(yk)))
           end do

        else if (idlink(yk).eq.2) then ! Splines link

           splaa=0.d0
           eta0=0.d0
           eta0=btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)

           do kk=2,ntr(yk)
              splaa(kk-3)=btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+kk)&
                   *btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+kk)
           end do
           !if(i==1 .and. id==0 .and. jd==0) print*,"eta0=",eta0,"splaa=",sqrt(splaa)
           do j=1,nmes(i,yk)
              ll=0
              !if(i==1 .and. id==0 .and. jd==0) print*,"Y=",Y(nmescur+sumMesYk+j)
              if (Y(nmescur+sumMesYk+j).eq.zitr(ntr(yk)-2,yk)) then
                 ll=ntr(yk)-3
              end if

              som=0.d0
              do kk = 2,ntr(yk)-2
                 if ((Y(nmescur+sumMesYk+j).ge.zitr(kk-1,yk)).and. &
                      (Y(nmescur+sumMesYk+j).lt.zitr(kk,yk))) then
                    ll=kk-1
                 end if
              end do

              if (ll.lt.1.or.ll.gt.ntr(yk)-3) then          
                 loglik=-1.d9
                 print*,"-1.d9 ll<1 ou ll>ntr-3",ll
                 goto 654
              end if
              if (ll.gt.1) then
                 do ii=2,ll
                    som=som+splaa(ii-3)
                 end do
              end if



              Y1(sumMesYk+j) = eta0 + som + &
                   splaa(ll-2) * im(2*ntotvalSPL+indiceY(nmescur+sumMesYk+j)) &
                   + splaa(ll-1) * im(ntotvalSPL+indiceY(nmescur+sumMesYk+j))&
                   + splaa(ll) * im(indiceY(nmescur+sumMesYk+j))

              jacobien = jacobien + log(splaa(ll-2)*mm(2*ntotvalSPL+indiceY(nmescur+sumMesYk+j)) &
                   +splaa(ll-1)*mm(ntotvalSPL+indiceY(nmescur+sumMesYk+j))&
                   +splaa(ll)*mm(indiceY(nmescur+sumMesYk+j)))
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
     !            write(*,*)'nmes',nmes(i),btot((nef+ncontr+nvc+1):npm),nef+ncontr
     !            write(*,*)'Y1',Y1
     !         end if




     ! contribution individuelle a la vraisemblance
     ! print*,"i=",i," -ni*log(2pi)=",-sum(nmes(i,:))*dlog(dble(2*3.14159265)), " log(det)=",det
     ! print*,"Vi=",VC
     vrais2_i=0.d0


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
              b0(l)=btot(nrisqtot+nvarxevt+nasso+l-1)
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
                 b01((m-1)*ny+yk)=btot(nrisqtot+nvarxevt+nasso+nef+(m-1)*(ny-1)+yk)
              else
                 b01((m-1)*ny+ny) =-sum(btot(nrisqtot+nvarxevt+nasso+nef+(m-1)*(ny-1)+1 &
                      :nef+(m-1)*(ny-1)+ny-1))
              end if
           end do
        end if
     end do


     som=0.d0
     som_T0=0.d0 !delayed entry
     som_Ti=0.d0
     do l=1,nMC
!print*,"Mc l=",l

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
        !print*,"mu=",mu

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
                    ai = abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
                 else if(methInteg.eq.2) then
                    !! MCA
                    if(mod(l,2).eq.0) then
                       ! si l est pair on prend l'oppose du precedent
                       ai = -abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
                    else
                       call bgos(SX,0,asim(yk),x22,0.d0)
                       ai = abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
                    end if
                 else
                    !! QMC
                    asim(yk) = seqMC(nMC*(nea+sum(nmes(i,:)))+l)
                    ai = abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk))*asim(yk)
                 end if
              end if
              !if(i.lt.4) print*,"i=",i," avant do j, vrais_Y=",vrais_Y
              do j=1,nmes(i,yk)

                 !! on ajoute ai a mu
                 if(nalea.gt.0) mu(sumMesYk+j) = mu(sumMesYk+j)+ai

                 !! trouver binf et bsup tq binf < lambda + epsilon < bsup

                 !! initialiser au premier seuil
                 binf = btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)
                 bsup = binf

                 !! si Y>minY ajouter b(..)^2
                 if(indiceY(nmescur+sumMesYk+j).gt.1) then
                    do ll=2,min(indiceY(nmescur+sumMesYk+j),ntr(yk))
                       bsup = bsup + btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+ll)**2
                       if(ll.lt.indiceY(nmescur+sumMesYk+j)) then
                          binf = binf + btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+ll)**2
                       end if
                    end do
                 end if
                 !if(i.lt.4)print*,"y=",Y1(sumMesYk+j)," indiceY=",indiceY(nmescur+sumMesYk+j), " Y=",Y(nmescur+sumMesYk+j)," mu=",mu(sumMesYk+j)
                 !print*," binf=",binf," bsup=",bsup
                 !! centrer et standardiser
                 binf = (binf - mu(sumMesYk+j))/abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk))
                 bsup = (bsup - mu(sumMesYk+j))/abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk))

                 !if(i.lt.4) print*,"j=",j," nvalORD=",nvalORD(ykord)," binf=",binf," bsup=",bsup

                 if(indiceY(nmescur+sumMesYk+j).eq.1) then
                    !! si Y=minY
                    vrais_Y = vrais_Y * alnorm(binf,.false.)
                    ! P(y=ymin) = P(y<=ymin) donc pareil si vrais ou expect
                 else if(indiceY(nmescur+sumMesYk+j).eq.nvalORD(ykord)) then
                    !! si Y=maxY

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
                    VC(j1,j1) = btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk)**2 !variance de l'erreur yk
                    if (nalea.eq.ny) then ! intercept aleatoire de yk
                       do j2=1,nmes(i,yk)
                          VC(j1,j2) = VC(j1,j2) + &
                               btot(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk)**2
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
                    loglik=-1.d9
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

!        if(i.lt.4) print*,"avant survie, vrais_Y",vrais_Y
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
                    brisq(k)=exp(btot(sumnrisq+k))
                 end do
              else
                 do k=1,nprisq(ke)
                    brisq(k)=btot(sumnrisq+k)*btot(sumnrisq+k)
                 end do
              end if

              basso=0.d0 
              basso=btot(nrisqtot+nvarxevt+ke) !basso = eta = prm estime des EAs partages

              beta_ef=0.d0 
              do k=1,nef
                 beta_ef(k) = btot(nrisqtot+nvarxevt+nasso+k)
              end do


              if (idst.eq.1) then  !bi
                 call fct_risq_loglik2(i,ke,brisq,risq,surv,surv0,survint, &
                      typrisq,logspecif,tsurv,zi,nz,idtrunc,tsurv0,nvdepsurv,tsurvint, &
                      Tmm,Tim,Tim0,Timt,ns,nbevt,nprisq)
              else if(idst.eq.2) then   !niv.courant du processus latent
!                 call fct_risq_irtsre_2(i,ke,brisq,basso,beta_ef,ui,risq,surv,surv0)  
              end if

              sumnrisq = sumnrisq + nprisq(ke)
           end do

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
                          bevt(m)=btot(nrisqtot+nevtxcurr+1)
                          Xevt(m)=X(nmescur+1,k)
                          nevtxcurr=nevtxcurr+1
                       else
                          if (idsurv(k).eq.2) then   
                             m=m+1
                             bevt(m)=btot(nrisqtot+nevtxcurr+ke)
                             Xevt(m)=X(nmescur+1,k)
                             nevtxcurr=nevtxcurr+nbevt
                          end if
                       end if

                    else ! i.e timedepvar

                       if (idsurv(k).eq.1) then  
                          bevtint(ke)=btot(nrisqtot+nevtxcurr+1)
                          nevtxcurr=nevtxcurr+1
                       else 
                          if (idsurv(k).eq.2) then
                             bevtint(ke)=btot(nrisqtot+nevtxcurr+ke)
                             nevtxcurr=nevtxcurr+nbevt
                          end if
                       end if
                    end if
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
                    easurv=DOT_PRODUCT(ui,btot((nrisqtot+nvarxevt+m+1)&
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
                    do ll=1,nef
                       pred_cl_Ti = pred_cl_Ti + Xcl_Ti(i,1+ll) * beta_EF(ll)  ! X(t) %*% beta
                    end do
                    do ll=1,nea
                       pred_cl_Ti = pred_cl_Ti + Xcl_Ti(i,1+nef+ll) * ui(ll)  ! Z(t) %*% ui
                    end do

                    pred_cl_Ti = EXP(pred_cl_Ti*btot(nrisqtot+nvarxevt+ke) )  ! exp(predcl*basso)
                    fevt=risq(ke)*exp(varexpsurv)*pred_cl_Ti
                 end if
                 if (nvdepsurv.eq.1) then
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

 !       if(l.lt.4) print*,"l=", l, " som=",som
        !    print*,"fin l=",l
     end do ! fin boucle nMC

  !print*,"fin MC vrais2_i=",vrais2_i, "jac=",jacobien
     vrais2_i = vrais2_i + log(som) - log(dble(nMC)) + jacobien

!     print*,"avt idtrunc vrais2_i=",vrais2_i
     if (idtrunc.eq.1) then    !delayedentry
        vrais2_i = vrais2_i - log(som_T0) + log(dble(nMC))
     end if


     loglik = loglik + vrais2_i

     if (vrais2_i.eq.-1.d9 .or. vrais2_i/vrais2_i.ne.1) then 
        loglik = -1.d9
        goto 654
     end if
     nmescur = nmescur + sum(nmes(i,:))
  end do

654 continue

  return

end subroutine loglik2








subroutine fct_risq_loglik2(i,k,brisq,risq,surv,surv0,survint, &
     typrisq,logspecif,tsurv,zi,nz,idtrunc,tsurv0,nvdepsurv,tsurvint, &
     Tmm,Tim,Tim0,Timt,ns,nbevt,nprisq)

        
  implicit none
  
  integer::i,k,ns,nbevt,logspecif,idtrunc,nvdepsurv
  double precision,dimension(nbevt)::risq,surv,surv0,survint
  double precision,dimension(ns)::tsurv,tsurv0,tsurvint
  integer,dimension(nbevt)::nz,typrisq,nprisq
  double precision,dimension(maxval(nz),nbevt)::zi
  double precision,dimension(nprisq(k))::brisq
  double precision,dimension(4*ns*nbevt)::Tmm,Tim,Tim0,Timt
  
  integer::j,l,ll,kk,ii
  double precision::som
  
  if (typrisq(k).eq.2.and.logspecif.eq.1) then
     
     surv(k)=brisq(1)*(tsurv(i)-zi(1,k))**brisq(2)      !zi(1,k)=depart Weibull
     
     risq(k)=brisq(1)*brisq(2)*(tsurv(i)-zi(1,k))**(brisq(2)-1)
     if (idtrunc.eq.1) then
        surv0(k)=brisq(1)*(tsurv0(i)-zi(1,k))**brisq(2)
     end if
     if (nvdepsurv.eq.1) then
        survint(k)=brisq(1)*(tsurvint(i)-zi(1,k))**brisq(2)
     else
        survint(k)=surv(k)
     end if
     
  end if
  if (typrisq(k).eq.2.and.logspecif.eq.0) then
     
     surv(k)=(brisq(1)*(tsurv(i)-zi(1,k)))**brisq(2)
     
     risq(k)=brisq(1)*brisq(2)*(brisq(1)*(tsurv(i)-zi(1,k)))**(brisq(2)-1)
     if (idtrunc.eq.1) then
        surv0(k)=(brisq(1)*(tsurv0(i)-zi(1,k)))**brisq(2)
     end if
     if (nvdepsurv.eq.1) then
        survint(k)=(brisq(1)*(tsurvint(i)-zi(1,k)))**brisq(2)
     else
        survint(k)=surv(k)
     end if
     
  end if
 
  if (typrisq(k).eq.1) then
     do j=1,nz(k)-1
        som=0.d0
        do l=1,j-1
           som=som+brisq(l)*(zi(l+1,k)-zi(l,k))
        end do
        if (idtrunc.eq.1) then
           if (Tsurv0(i).ge.zi(j,k).and.Tsurv0(i).le.zi(j+1,k)) then
              surv0(k)=som+brisq(j)*(Tsurv0(i)-zi(j,k))
           end if
        end if
        if (Tsurv(i).ge.zi(j,k).and.Tsurv(i).le.zi(j+1,k)) then
           surv(k)=som+brisq(j)*(Tsurv(i)-zi(j,k))
           risq(k)=brisq(j)
        end if
        if (nvdepsurv.eq.1) then
           if (Tsurvint(i).ge.zi(j,k).and.Tsurvint(i).le.zi(j+1,k)) &
                then
              survint(k)=som+brisq(j)*(Tsurvint(i)-zi(j,k))
           end if
        end if
     end do
     if (nvdepsurv.eq.0) then
        survint(k)=surv(k)
     end if
  end if
  
  

  if (typrisq(k).eq.3) then
     !------------ survie et risq pour Tsurv ----------------
     ll=0
     if (Tsurv(i).eq.zi(nz(k),k)) then
        ll=nz(k)-1
     end if
     som=0.d0
     do kk=2,nz(k)
        if ((Tsurv(i).ge.zi(kk-1,k)).and.(Tsurv(i).lt.zi(kk,k))) &
             then
           ll=kk-1
        end if
     end do
     if (ll.gt.1) then
        do ii=1,ll-1
           som=som+brisq(ii)
        end do
     end if
     
     surv(k) = som + brisq(ll) * Tim(3*ns*nbevt+ns*(k-1)+i) + &
          brisq(ll+1) * Tim(2*ns*nbevt+ns*(k-1)+i) + &
          brisq(ll+2) * Tim(ns*nbevt+ns*(k-1)+i) + &
          brisq(ll+3) * Tim(ns*(k-1)+i)
     risq(k) = brisq(ll) * Tmm(3*ns*nbevt+ns*(k-1)+i) + &
          brisq(ll+1) * Tmm(2*ns*nbevt+ns*(k-1)+i) + &
          brisq(ll+2) * Tmm(ns*nbevt+ns*(k-1)+i) + &
          brisq(ll+3) * Tmm(ns*(k-1)+i)

     !------------ survie et risq pour Tsurv0 ----------------
     
     if (idtrunc.eq.1) then
        ll=0
        if (Tsurv0(i).eq.zi(nz(k),k)) then
           ll=nz(k)-1
        end if
        som=0.d0
        do kk=2,nz(k)
           if ((Tsurv0(i).ge.zi(kk-1,k)).and.(Tsurv0(i).lt.zi(kk,k))) &
                then
              ll=kk-1
           end if
        end do
        !               if (ll.lt.1.or.ll.gt.nz-1) then
        !                  write(*,*) 'probleme dans fct_risq splines'
        !                  write(*,*) 'll=',ll,'T=',Tsurv0(i)
        !                  stop
        !               end if
        if (ll.gt.1) then
           do ii=1,ll-1
              som=som+brisq(ii)
           end do
        end if
        
        surv0(k) = som + brisq(ll) * Tim0(3*ns*nbevt+ns*(k-1)+i) + &
             brisq(ll+1) * Tim0(2*ns*nbevt+ns*(k-1)+i) + &
             brisq(ll+2) * Tim0(ns*nbevt+ns*(k-1)+i) + &
             brisq(ll+3) * Tim0(ns*(k-1)+i)
        
     end if
     
     !------------ survie et risq pour Tsurvint ----------------


     if (nvdepsurv.eq.1) then

        !               write(*,*)'i',i,tsurvint(i),ind_survint(i),tsurv(i),tsurv0(i)
!               write(*,*)timt3(i),Timt2(i),timt1(i)

        ll=0
        if (Tsurvint(i).eq.zi(nz(k),k)) then
           ll=nz(k)-1
        end if
        som=0.d0
        do kk=2,nz(k)
           if((Tsurvint(i).ge.zi(kk-1,k)).and.(Tsurvint(i).lt.zi(kk,k))) &
                then
              ll=kk-1
           end if
        end do
        !               if (ll.lt.1.or.ll.gt.nz-1) then
        !                  write(*,*) 'probleme dans fct_risq splines'
        !                  write(*,*) 'll=',ll,'T=',Tsurvint(i)
        !                  stop
        !               end if
        if (ll.gt.1) then
           do ii=1,ll-1
              som=som+brisq(ii)
           end do
        end if
               
        survint(k) = som + brisq(ll) * Timt(3*ns*nbevt+ns*(k-1)+i)+ &
             brisq(ll+1) * Timt(2*ns*nbevt+ns*(k-1)+i) + &
             brisq(ll+2) * Timt(ns*nbevt+ns*(k-1)+i) + &
             brisq(ll+3) * Timt(ns*(k-1)+i)
        
     else
        survint(k) = surv(k)
     end if
     
  end if
  
  
end subroutine fct_risq_loglik2
