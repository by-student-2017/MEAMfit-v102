subroutine amendFreep

    !---------------------------------------------------------------c
    !
    !     This subroutine checks for non-zero elements of the 
    !     potential file being read in, and if it finds them, and
    !     if the corresponding freep value is set to 2 (initialize
    !     parameter from random seed), then it will reset freep to
    !     0 or 1, so that the value from the file is instead read 
    !     in.
    !
    !     Returns:       freep
    !
    !     Andrew Duff
    !
    !---------------------------------------------------------------c

    use m_optimization
    use m_atomproperties

    implicit none

    integer i

    do i=1,m4+2*m2+m1
        if ((p(i).ne.0d0).and.(freep(i).eq.2)) then
!           print *,'Changing freep(',i,')=2 to freep(',i,')=1'
           if (fixPotIn.eqv..true.) then
              freep(i)=0
           else
              freep(i)=1
           endif
        endif
    enddo

end subroutine amendFreep
subroutine backgrounddensity

    !--------------------------------------------------------------c
    !
    !     Calculates the background densities rho_i for each of
    !     the inequivalent sites (these are then used to calculate
    !     the embedding functions).
    !
    !     Called by:     meamenergy
    !     Calls:         -
    !     Arguments:     istr,gn_inequivalentsites,meam_t,rhol
    !     Returns:       rho_i
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter and Andrew Duff Feb 9 2006
    !
    !--------------------------------------------------------------c

    use m_generalinfo
    use m_geometry
    use m_meamparameters
    use m_electrondensity

    implicit none

    integer i,angmom
    real*8 gamma

    if(allocated(rho_i)) deallocate(rho_i)
    allocate(rho_i(gn_inequivalentsites(istr)))

    !Loop over all atom sites for atomic configuration, istr
    do i=1,gn_inequivalentsites(istr)
        !Construct gamma, to be used in the equation for rho_i
        gamma=0d0
        do angmom=1,lmax
            if (envdepmeamt.eqv..false.) then
                !Angular dependant densities equal directly to the
                !meamtau values specified in the meam parameters file
                !(rather than the averaged quantity which is used in
                !Camelion, for example)
                meam_t(angmom,i)=meamtau(angmom,gspecies(i,istr))
            endif
            gamma=gamma+meam_t(angmom,i)*(rhol(angmom,i)**2)
        enddo
        if (rhol(0,i).ne.0d0) then
            gamma=gamma/(rhol(0,i)**2)
        endif

        !Calculate total electron density on site i
        rho_i(i)=rhol(0,i)*2d0/(1d0+exp(-gamma))

    enddo

end subroutine backgrounddensity
subroutine backupxyz

    !Backup xyz files. This will be used once we have structural relaxation
    !included: After each optimization the positions can then be reset to their
    !original positions
    !
    !Andrew Duff, 2013

    use m_geometry

    integer i,j

    !Store initial atom positions in gxyz_backup
    do i=1,nstruct
        do j=1,gn_inequivalentsites(i)
            gxyz_backup(1:3,j,i)=gxyz(1:3,j,i)
        enddo
    enddo

end subroutine backupxyz
subroutine broyden_bns(iter,mixparam,nbro,nstor, &
        q,qinv,f_old,aux,aux2,y,s, &
        x,f)
    ! multi-secant Broyden method to solve f(x)=0
    ! storage saving BNS algorithm
    ! first iteration is linear mixing with mixparam
    ! multi-secant condition gives singular matrix Q, therefore SVD is used
    ! algorithm is modification of Sawamura et al (Trans JIM 40, p1186,
    ! 1999)
    ! Marcel Sluiter, April 24 2000
    implicit none
    integer iter,nbro,nstor,istor,istorp,istorn,i,j,idim
    real(8) x(nbro),f(nbro),mixparam,nrmy, &
        s(nbro,nstor),y(nbro,nstor),q(nstor,nstor),qinv(nstor,nstor), &
        f_old(nbro),aux(nstor),aux2(nstor)

    istorp=mod(iter-3,nstor)+1   !previous storage location
    istor =mod(iter-2,nstor)+1   !current storage location
    istorn=mod(iter-1,nstor)+1   !next storage location
    if(iter.ge.2) then
        y(1:nbro,istor) = f - f_old
        nrmy = sqrt(sum(y(1:nbro,istor)**2))
        y(1:nbro,istor) = y(1:nbro,istor) / nrmy
        s(1:nbro,istor) = s(1:nbro,istor) / nrmy
    endif
    f_old = f

    ! generate Q-matrix
    idim = min(iter-1,nstor)
    do j = 1, idim
        q(istor,j) = -dot_product(y(:,istor),y(:,j))
        q(j,istor) = q(istor,j)
    enddo
    if(idim.ge.2) then
        ! get q**(-1) using SVD for idim=>2
        call svdinv(q,idim,idim,nstor,nstor, &
            qinv)
    else if(idim.eq.1) then
        if(abs(q(1,1)).gt.1d-6) then
            qinv(1,1)=1d0/q(1,1)
        else
            qinv(1,1)=1d0
        endif
    endif

    ! compute x(n+1) = x(n) + b[I+y(n-1)(q(n-1)**-1)y(n-1)]f(n) +
    !                       + s(n-1)(q(n-1)**-1)y(n-1)f(n)
    ! where y(j) = f(j+1)-f(j), s(j) = x(j+1)-x(j)
    ! aux = y(n-1)*f(n)
    do i = 1, idim
        aux(i) = dot_product(y(1:nbro,i),f)
    enddo
    ! aux2 = (q(n-1)**-1)y(n-1)f(n) = (q(n-1)**-1) aux
    do i = 1, idim
        aux2(i) = dot_product(qinv(i,1:idim),aux(1:idim))
    enddo
    ! xnew temporarily stored in f
    ! xnew = f(n) + y(n-1)(q(n-1)**-1)y(n-1)f(n) = f(n) + y(n-1) aux2
    do i = 1, idim
        f(1:nbro) = f(1:nbro) + y(1:nbro,i)*aux2(i)
    enddo
    ! xnew = b * [f(n) + y(n-1)(q(n-1)**-1)y(n-1)f(n) ]
    f = mixparam * f
    ! xnew = b * [f(n) + y(n-1)(q(n-1)**-1)y(n-1)f(n) ] +
    !                    s(n-1)(q(n-1)**-1)y(n-1)f(n)
    do i = 1, idim
        f(1:nbro) = f(1:nbro) + s(1:nbro,i)*aux2(i)
    enddo
    ! x(n+1) = x(n) + b (I + y(n-1).....
    x = x + f
    s(1:nbro,istorn) = f
end subroutine broyden_bns
subroutine broyden(iter,nbro,nstor,mixparam, &
        work, &
        x,f)
    ! driver program for Broyden, solves f(x)=0
    ! give a trial x and f(x), returns new guess for x
    ! needs to be called repeatedly until ||f|| is sufficiently small
    ! Marcel Sluiter, April 26 2000, Nov 28 2001

    ! iter  = current iteration number (first call, iter=1)
    ! nbro  = dimensionality of problem
    ! nstor = number of iterations that are stored
    ! mixparam = mixing parameter for first step. [0:1], usually about 0.1
    ! work     = work area: real(8) dimension(nbro+2*nstor*(1+nstor+nbro))
    !            DO NOT overwrite work!
    ! x        = input: guess for solving f(x)=0
    !            output: new guess
    ! f        = input: function to be set to zero, evaluated at x
    !            output: x(new_guess) - x(old_guess), change in x

    implicit none
    integer iter,nbro,nstor,n,m
    real(8) x(nbro),f(nbro),mixparam,work(nbro+2*nstor*(1+nstor+nbro))
    intent(in) iter,nbro,nstor,mixparam

    n=nstor**2
    m=1+2*n+nbro+2*nstor
    call broyden_bns(iter,mixparam,nbro,nstor, &
        work(1),work(1+n),work(1+2*n),work(1+2*n+nbro), &
        work(1+2*n+nbro+nstor),work(m),work(m+nbro*nstor), &
        x,f)
end subroutine broyden
subroutine CalcCutoffPenalty

    !----------------------------------------------------------------------c
    !
    !     Currently un-used subroutine, but may be re-implemented again at 
    !     later date. Designed to apply a penalty to the optimization 
    !     function proportional to the square of the deviation of the 
    !     potential parameters from their starting values (originally to 
    !     try to control the cutoff radii, to stop them moving 
    !     uncontrollably beyond their limiting values).
    !
    !     Called by:     originally by 'fit_quality'
    !     Returns:       CutoffPenalty
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Oct 2014
    !
    !----------------------------------------------------------------------c

    use m_optimization
    use m_generalinfo
    use m_meamparameters
    use m_atomproperties

    implicit none

    integer i,ip,ispc,jspc,spc_cnt,ll

    CutoffPenalty=0d0
    CutoffPenCoeff=1d0 !Determines severity of penalty for deviation

    !Get contributions to CutoffPenalty from electron density cutoffs
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    do i=2,12,2
                       ip=2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i
                       CutoffPenalty=CutoffPenalty+CutoffPenCoeff*((p(ip)-p_orig(ip))**2)
                    enddo
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo
    
    !!Get contributions to CutoffPenalty from pair-potential cutoffs
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                do i=2,32,2
                   ip=2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt
                   CutoffPenalty=CutoffPenalty+CutoffPenCoeff*((p(ip)-p_orig(ip))**2)
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

end subroutine CalcCutoffPenalty
subroutine calcSds

    !----------------------------------------------------------------------c
    !
    !     Calculate the standard-deviations of the input energies and
    !     forces for rescaling of the optimization function.
    !
    !     Called by:     program MEAMfit
    !     Returns:       
    !     Files read:    
    !     Files written: 
    !
    !     Andy Duff, Oct 2014
    !
    !----------------------------------------------------------------------c

    use m_geometry
    use m_datapoints
    use m_optimization

    implicit none

    integer i,istruc,idatapoint

    !Calculate averages of energies and forces
    idatapoint=1
    avgEn=0d0
    avgFrc=0d0
    nEn=0
    nFrcComp=0

    do istruc=1,nstruct
        if (optimizeforce(istruc).gt.0) then
            do i=1,gn_forces(istruc)
               if (optforce(i,istruc).eqv..true.) then
                  avgFrc = avgFrc + truedata(idatapoint) + truedata(idatapoint+1) + truedata(idatapoint+2)
                  nFrcComp=nFrcComp+3
               endif
               idatapoint=idatapoint+3
            enddo
        else
            avgEn=avgEn+truedata(idatapoint)
            idatapoint=idatapoint+1
            nEn=nEn+1
        endif
    enddo
    if (nEn.gt.0) then
       avgEn=avgEn/nEn
    endif
    if (nFrcComp.gt.0) then
       avgFrc=avgFrc/nFrcComp
    endif

    !Calculate variances of energies and forces
    idatapoint=1
    varEn=0d0
    varFrc=0d0
    do istruc=1,nstruct
        if (optimizeforce(istruc).gt.0) then
            do i=1,gn_forces(istruc)
               if (optforce(i,istruc).eqv..true.) then
                  varFrc=varFrc+(truedata(idatapoint)-avgFrc)**2 + &
                                (truedata(idatapoint+1)-avgFrc)**2 + &
                                (truedata(idatapoint+2)-avgFrc)**2
         !     print *,'idatapoint=',idatapoint,', truedata=',truedata(idatapoint),' diff=', &
         !             (truedata(idatapoint)-avgFrc)**2
               endif
               idatapoint=idatapoint+3
            enddo
        else
            varEn=varEn+(truedata(idatapoint)-avgEn)**2
            idatapoint=idatapoint+1
        endif
    enddo
    if (nEn.gt.0) then
       varEn=varEn/nEn
    endif
    if (nFrcComp.gt.0) then
       varFrc=varFrc/nFrcComp
       !print *,'varFrc before division=',varFrc,', nFrcComp=',nFrcComp
    endif

    !print *,'idatapoint=',idatapoint,', ndatapoints=',ndatapoints
    print *,'avgEn=',avgEn,', avgFrc=',avgFrc
    print *,'varEn=',varEn,', varFrc=',varFrc

end subroutine calcSds
subroutine calculateshells(nshell,nnatoms_est,rmax)

    !----------------------------------------------------------------------c
    !
    !     Calculates how many shells, 'nshell', of unit cells are required
    !     in order to include all atoms which are within a distance rmax
    !     from at least one of the atoms in the
    !     central unit cell. The (exact) number of such atoms
    !     atoms is calculated as 'nnatoms_est'.
    !     The shells are defined so that, for example, the first shell of
    !     unit cells forms a cubic shell in direct coordinate space about
    !     the central unit cell. I.e., the unit cells included in the
    !     first shell are those at (in direct coordinates): (-1,-1,-1),
    !     (0,-1,-1), (1,-1,-1), (-1,0,-1), (0,0,-1), (1,0,-1), (-1,1,-1),
    !     (0,1,-1), (1,1,-1) <- bottom face; (-1,-1,0), (0,-1,0),
    !     (1,-1,0), (-1,0,0), (1,0,0), (-1,1,0), (0,1,0), (1,1,0) -< mid-
    !     section; (-1,-1,1), (0,-1,1), (1,-1,1), (-1,0,1), (0,0,1),
    !     (1,0,1), (-1,1,1) -<top face.
    !
    !     Called by:     neighborhood
    !     Calls:         map_icoords_to_its
    !     Returns:       nshell,nnatoms_est
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Dec 2007
    !
    !----------------------------------------------------------------------c

    use m_generalinfo
    use m_poscar
    use m_neighborlist

    implicit none
    logical no_nn_found,nn_found
    logical, allocatable:: atom_tested(:,:,:,:) !Use this array to
    !record which atoms have been tested to see if they
    !are nn's. This is so that we don't triple count
    !the number of nn's in the cells of adjoining
    !faces. The first 3 indicies are for the unit cell,
    !the fourth is for the atom in the unit cell.
    integer j,k,it1,it2,it3,iface,icoord1,icoord2,nshell,nnatoms_est
    integer, parameter:: unitcell_dim=10 !The maximum number of unit
    !cells expected out from the origin along any
    !lattice direction.
    real(8) rmax2,xyztry(3),dxyz(3),dd,rmax

    allocate(atom_tested(-unitcell_dim:unitcell_dim, &
        -unitcell_dim:unitcell_dim,-unitcell_dim:unitcell_dim, &
        natoms))

    atom_tested=.false.
    rmax2=rmax**2
    nshell=1
    nnatoms_est=n_inequivalentsites

    do                        !Loop over 'shells' of unit cells around
        no_nn_found=.true.     !the central unit cell

        do iface=1,6           !Loop over the faces of the cube (a cube
            !in direct coordinate space)

            do icoord1=-nshell,nshell !Loop over the points on the
                do icoord2=-nshell,nshell !cube face

                    call map_icoords_to_its(it1,it2,it3,iface,nshell, &
                        icoord1,icoord2) !Convert icoord values into
                    !it values.

                    do j=1,natoms

                        xyztry=coordinates(1:3,j)+tr(1:3,1)*real(it1)+ &
                            tr(1:3,2)*real(it2)+tr(1:3,3)*real(it3)

                        !Check to see if this coordinate is a nearest
                        !neighbor of any of the atoms in the central unit
                        !cell
                        nn_found=.false.
                        do k=1,n_inequivalentsites
                            dxyz=xyztry-coordinates(1:3,k)
                            dd=dxyz(1)**2+dxyz(2)**2+dxyz(3)**2

                            if (dd.le.rmax2) then
                                nn_found=.true.
                                no_nn_found=.false.
                            endif

                        enddo

                        !The same unit cells are included more than once
                        !on the separate faces. Take this into account in
                        !the estimate for the number of nn's:
                        if ((nn_found.eqv..true.).and. &
                            (atom_tested(it1,it2,it3,j).eqv..false.)) &
                            then
                            nnatoms_est=nnatoms_est+1
                        endif
                        atom_tested(it1,it2,it3,j)=.true.

                    enddo
                    !
                enddo
            enddo
        enddo

        if (no_nn_found) then  !No nn's found the present shell,
            exit                !Therefore we have located all nn's of
        endif                  !the central unit cell
        nshell=nshell+1
    enddo

    nshell=nshell-1
    deallocate(atom_tested)

end subroutine calculateshells
subroutine calculateshells_forces(nshell,nnatoms_est,rmax)

    !----------------------------------------------------------------------c
    !
    !     Calculates how many shells, 'nshell', of unit cells are required
    !     in order to include all atoms which are within a distance rmax
    !     from at least one of the atoms in the
    !     central unit cell. The (exact) number of such atoms
    !     atoms is calculated as 'nnatoms_est'.
    !     The shells are defined so that, for example, the first shell of
    !     unit cells forms a cubic shell in direct coordinate space about
    !     the central unit cell. I.e., the unit cells included in the
    !     first shell are those at (in direct coordinates): (-1,-1,-1),
    !     (0,-1,-1), (1,-1,-1), (-1,0,-1), (0,0,-1), (1,0,-1), (-1,1,-1),
    !     (0,1,-1), (1,1,-1) <- bottom face; (-1,-1,0), (0,-1,0),
    !     (1,-1,0), (-1,0,0), (1,0,0), (-1,1,0), (0,1,0), (1,1,0) -< mid-
    !     section; (-1,-1,1), (0,-1,1), (1,-1,1), (-1,0,1), (0,0,1),
    !     (1,0,1), (-1,1,1), (0,1,1), (1,1,1) -<top face.
    !
    !     Note: more atoms are included than in 'calculateshells', since
    !     an atom is displaced in the central unit cell
    !
    !     Called by:     neighborhood
    !     Calls:         map_icoords_to_its
    !     Returns:       nshell,nnatoms_est
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Dec 2007
    !
    !----------------------------------------------------------------------c

    use m_generalinfo
    use m_poscar
    use m_neighborlist

    implicit none
    logical inlist,no_nn_found,nn_found
    logical, allocatable:: atom_tested(:,:,:,:) !Use this array to
    !record which atoms have been tested to see if they
    !are nn's. This is so that we don't triple count
    !the number of nn's in the cells of adjoining
    !faces. The first 3 indicies are for the unit cell,
    !the fourth is for the atom in the unit cell.
    integer j,k,it1,it2,it3,iface,icoord1,icoord2,nshell, &
        nnatoms_est,atomnumber
    integer, parameter:: unitcell_dim=10 !The maximum number of unit
    !cells expected out from the origin along any
    !lattice direction.
    real(8) rmax2,xyztry(3),dxyz(3),dd,rmax

    allocate(atom_tested(-unitcell_dim:unitcell_dim, &
        -unitcell_dim:unitcell_dim,-unitcell_dim:unitcell_dim, &
        natoms))

    atom_tested=.false.
    rmax2=rmax**2
    nshell=1
    nnatoms_est=n_inequivalentsites

    do                        !Loop over 'shells' of unit cells around
        no_nn_found=.true.     !the central unit cell

        do iface=1,6           !Loop over the faces of the cube (a cube
            !in direct coordinate space)

            do icoord1=-nshell,nshell !Loop over the points on the
                do icoord2=-nshell,nshell !cube face

                    call map_icoords_to_its(it1,it2,it3,iface,nshell, &
                        icoord1,icoord2) !Convert icoord values into
                    !it values.

                    do j=1,natoms_old

                        xyztry=coords_old(1:3,j)+tr(1:3,1)*real(it1)+ &
                            tr(1:3,2)*real(it2)+tr(1:3,3)*real(it3)

                        !Check to see if this coordinate is a nearest
                        !neighbor of any of the inequivalent atoms
                        nn_found=.false.
                        do k=1,n_inequivalentsites
                            dxyz=xyztry-coordinates(1:3,k)
                            dd=dxyz(1)**2+dxyz(2)**2+dxyz(3)**2

                            if (dd.le.rmax2) then
                                !First check that the nearest neighbor isn't
                                !actually an inequivalent site itself
                                call checkatom(xyztry,inlist,atomnumber)
                                if (inlist.eqv..false.) then
                                    nn_found=.true.
                                    no_nn_found=.false.
                                endif
                            endif

                        enddo

                        !The same unit cells are included more than once
                        !on the separate faces. Take this into account in
                        !the estimate for the number of nn's:
                        if ((nn_found.eqv..true.).and. &
                            (atom_tested(it1,it2,it3,j).eqv..false.)) &
                            then
                            nnatoms_est=nnatoms_est+1
                        endif
                        atom_tested(it1,it2,it3,j)=.true.

                    enddo
                    !
                enddo
            enddo
        enddo

        if (no_nn_found) then  !No nn's found the present shell,
            exit                !Therefore we have located all nn's of
        endif                  !the central unit cell
        nshell=nshell+1
    enddo

    nshell=nshell-1
    if (nshell.eq.0) then
        print *,'calculateshells_forces should never return'
        print *,'nshell=0'
    endif
    deallocate(atom_tested)

end subroutine calculateshells_forces
subroutine checkatom(xyztry,inlist,atomnumber)

    !--------------------------------------------------------------c
    !
    !     Takes xyztry and finds out if it is equal to xyz(1:3,j) 
    !     of any of the inequivalent atoms. If so, inlist is set to 
    !     .true. and the atomnumber of the atom is returned
    !
    !     Called by:     neighborhood
    !     Calls:         -
    !     Returns:       inlist, atomnumber
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff, 2010
    !
    !--------------------------------------------------------------c

    use m_poscar

    implicit none

    integer atomnumber,i
    real(8) xyztry(3)
    logical inlist

    inlist=.false.
    do i=1,natoms
        !Check components of xyztry against those of inequivalent atoms
        if ( (coordinates(1,i).eq.xyztry(1)).and. &
            (coordinates(2,i).eq.xyztry(2)).and. &
            (coordinates(3,i).eq.xyztry(3)) ) then
            atomnumber=i
            inlist=.true.
            exit
        endif
    enddo

end subroutine checkatom
subroutine checkInputFilesExist

    !--------------------------------------------------------------c
    !
    !     Check to see if a settings file is present. Also, if the
    !     user has specified for a potential file to be read in,
    !     check it exists.
    !
    !     Called by:     program MEAMfit
    !     Returns:       settingsfileexist
    !
    !     Andrew Duff 2014
    !
    !--------------------------------------------------------------c

    use m_filenames
    use m_generalinfo

    implicit none

    logical exist

    !If a potential parameter file has been specified to be read-in, make sure
    !it exists.
    if (readpotfile.eqv..true.) then
       inquire(file=trim(startparameterfile),exist=exist)
       if (exist.eqv..false.) then
          print *,'ERROR: specified read in potential parameters file ', &
                 trim(startparameterfile),' but file does not exist, stopping.'
          stop
       endif
    endif

    !Check if settings file exists
    inquire(file=trim(settingsfile),exist=exist)
    if (exist.eqv..false.) then
       settingsfileexist=.false.
    else
       settingsfileexist=.true.
    endif


end subroutine checkInputFilesExist
subroutine checkParas

    !----------------------------------------------------------------------c
    !
    !     Check potential parameters and variables used to set up these
    !     parameters. Display information about potential to user.
    !
    !     Called by: MEAMfit
    !     Calls: -
    !     Returns: cmin_cmax_zero, pairpotonly 
    !     Files read:   
    !     Files written:
    !
    !     Andrew Duff, 2014
    !
    !----------------------------------------------------------------------c

    use m_meamparameters
    use m_generalinfo
    use m_atomproperties
    use m_optimization
    use m_filenames

    implicit none

    logical cutoffsOk,firstoccurance
    integer i,j,k,ispc,jspc,ll,spc_cnt

    !Check to see if all the cmin and cmax equal zero and if they are
    !not free parameters.
    cmin_cmax_zero=.true.
    if ((freep(1).eq.1).or.(freep(2).eq.1).or.(freep(1).eq.2).or.(freep(2).eq.2)) then
        cmin_cmax_zero=.false.
    endif
    do i=1,maxspecies
        do j=1,maxspecies
            do k=1,maxspecies
                if ((cmin(i,j,k).ne.0).and.(cmax(i,j,k).ne.0)) then
                    cmin_cmax_zero=.false.
                endif
            enddo
        enddo
    enddo
    if (cmin_cmax_zero.eqv..false.) then
        print *,'Baskes angular screening is on'
    endif

    !Check if we have only a pair-potential
    if ((cmin_cmax_zero.eqv..true.).and.(typethi(1,1).eq.0).and. &
        (typethi(1,2).eq.0).and.(typethi(2,2).eq.0).and. &
        (embfunctype.eq.0)) then
        pairpotonly=.true.
        print *,'Pair-potential only'
        print *,'WARNING: output file for LAMMPS needs to be adjusted: stopping)'
        stop
    endif

    if (cutoffopt.eqv..true.) then
       print *,'Cut-off radii are being optimized, within range CUTOFF_MIN ', &
              '< cut-off radius < CUTOFF_MAX'
       print *
       print *,'Checking cut-off radii:'
       !Check that p_rmax is larger than all of the cut-offs which are to be
       !initialized from the start_meam_parameters file. Also, for those cut-offs 
       !which are to be randomly generated, check that p_rmax is larger than the 
       !corresponding upper limits for the random generation.
       cutoffsOk=.true.
       spc_cnt=0
       do ispc=1,maxspecies
           do jspc=1,maxspecies
               if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                   do ll=0,lmax
                       !Check upper-limits for the cut-off random generation
                       !(those upper-limits not in use were set to zero in 
                       !initBoundsForRandParas, so will not wrongly trigger 
                       !the following if/then condition)
                       if (meamrhodecay_maxradius(ll,ispc,jspc).gt.p_rmax) then
                           if (cutoffsOk.eqv..true.) then
                              print *,'ERROR:'
                           endif
                           print *,'The upper-limit for random-generation of the electron density cut-off'
                           write(*,'(A17,I1,A7,I1,A4,I1,A24)') ' radii (for ispc=',ispc,', jspc=',jspc, &
                                      ', l=',ll,') is larger than p_rmax.'
                           cutoffsOk=.false.
                       endif
                       !Check cut-offs which are not to be randomly generated
                       do i=2,12,2
                           if ((freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i) &
                                .lt.2).and.(meamrhodecay(i,ll,ispc,jspc) &
                                .gt.p_rmax)) then
                               if (cutoffsOk.eqv..true.) then
                                  print *,'ERROR:'
                               endif
                               write(*,'(A46,A10,I1,A17)') ' The starting value of the cut-off radius for ', & 
                                    'parameter ',i,' of the electron-'
                               write(*,'(A18,I1,A7,I1,A4,I1,A24)') ' density (for ispc=',ispc, &
                                    ', jspc=',jspc,', l=',ll,') is larger than p_rmax.'
                               cutoffsOk=.false.
                           endif
                       enddo
                  enddo
              endif
              spc_cnt=spc_cnt+1
           enddo
        enddo

       spc_cnt=0
       do ispc=1,maxspecies
           do jspc=1,maxspecies
               if (jspc.ge.ispc) then
                   !Check that upper-limit for cut-off random generation is
                   !less than p_rmax
                   if (pairpotparameter_maxradius(ispc,jspc).gt. &
                       p_rmax) then
                       if (cutoffsOk.eqv..true.) then
                          print *,'ERROR:'
                       endif
                       print *,'The upper-limit for random-generation of the pair-potential cut-off'
                       write(*,'(A17,I1,A7,I1,A24)') ' radii (for ispc=',ispc,', jspc=',jspc, &
                                  ') is larger than p_rmax.'
                       cutoffsOk=.false.
                   endif
                   !Also check that those cut-offs which are not to be randomly
                   !generated are also less than p_rmax
                   do i=2,32,2
                       if ((freep(2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt) &
                           .lt.2).and.(pairpotparameter(i,ispc,jspc).gt.p_rmax)) then
                           if (cutoffsOk.eqv..true.) then
                              print *,'ERROR:'
                           endif
                           write(*,'(A46,A10,I1,A13)') ' The starting value of the cut-off radius for ', &
                                'parameter ',i,' of the pair-'
                           write(*,'(A20,I1,A7,I1,A24)') ' potential (for ispc=',ispc, &
                                ', jspc=',jspc,') is larger than p_rmax.'
                           cutoffsOk=.false.
                       endif
                   enddo
               endif
               spc_cnt=spc_cnt+1
           enddo
        enddo

        if (cutoffsOk.eqv..false.) then
           print *
           print *,'Please increase p_rmax or decrease the values of the cutoffs'
           print *,'(or their upper limits) and run again.'
           print *
           print *,'STOPPING.'
           stop
        else
           print *,'...all good'
        endif

    endif
    !--------------------------------------------------------------------


    !---- Check if any finite value potential-parameters have not been set to optimize ----
    !If so: warn the user. Note that these warnings are less relevant than when
    !I originally wrote the code, and so are currently only produced when
    !verbose=.true.
    if ((nsteps.gt.1).and.verbose) then

        firstoccurance=.true.

        !Check the coefficients of the electron-density
        do i=2*m3+lm1*m1+1,2*m3+lm1*m1+12*lm1*m2,2
           !If acceptor independant densities, then only read half this array
           if (i.le.2*m3+lm1*m1+12*lm1*m1) then
              if ((freep(i).eq.0).and.(p(i).ne.0d0).and.(p(i+1).ne.0d0)) then
                 if (firstoccurance) then
                    print *
                    print *,'WARNING:'
                    firstoccurance=.false.
                 endif
                 write(*,*) 'You have a coefficient/cutoff radius pair in your electron-density (parameter'
                 write(*,'(A1,I2,A5,I2,A53)') ' ',i,' and ',i+1,' in the start_meam_parameters file) which is non-zero'
                 write(*,*) '(',p(i),',',p(i+1),')'
                 write(*,'(A61)') 'and yet which is set to remain fixed during the optimization'
              endif
           endif
        enddo

        !Check the coefficients of the embedding function
        do i=2*m3+lm1*m1+12*lm1*m2+1,2*m3+4*m1+lm1*m1+12*lm1*m2
           if ((freep(i).eq.0).and.(p(i).ne.0d0)) then
               if (firstoccurance) then
                  print *
                  print *,'WARNING:'
                  firstoccurance=.false.
               endif
               write(*,'(A53,I3,A8)') ' You have an embedding function parameter (parameter ',i,') in the'
               write(*,'(A46)') ' start_meam_parameters file) which is non-zero'
               write(*,*) '(',p(i),')'
               write(*,'(A61)') 'and yet which is set to remain fixed during the optimization'
            endif
        enddo 

        !Check the coefficients of the pair-potential
        spc_cnt=0
        do ispc=1,maxspecies
            do jspc=1,maxspecies
                if (jspc.ge.ispc) then
                   do i=2*m3+(4+lm1)*m1+12*lm1*m2+1+spc_cnt*32,2*m3+(4+lm1)*m1+12*lm1*m2+32+spc_cnt*32,2
                      if ((freep(i).eq.0).and.(p(i).ne.0d0).and.(p(i+1).ne.0d0)) then
                         if (firstoccurance) then
                            print *
                            print *,'WARNING:'
                            firstoccurance=.false.
                         endif
                         write(*,*) 'You have a coefficient/cutoff radius pair in your pair-potential (parameter'
                         write(*,'(A1,I3,A5,I3,A53)') ' ',i,' and ',i+1,' in the start_meam_parameters file) which is non-zero'
                         write(*,*) '(',p(i),',',p(i+1),')'
                         write(*,'(A61)') 'and yet which is set to remain fixed during the optimization'
                      endif
                   enddo
                endif
                spc_cnt=spc_cnt+1
            enddo
        enddo

        if (firstoccurance.eqv..false.) then
           print *
           print *,'...is this what you intended?'
        endif

    endif
    !--------------------------------------------------------------------------------------

end subroutine checkParas
subroutine createfitdbse

    !----------------------------------------------------------------------c
    !
    !     Checks for vasprun.xml files in directory (later add in OUTCAR
    !     file check) and makes a template fitdbse file from them.
    !
    !     Called by: setupfilenames    
    !     Calls: system
    !     Files read: - 
    !     Files written: fitdbse
    !
    !     Andrew Duff, 2014
    !
    !----------------------------------------------------------------------c

    use m_meamparameters
    use m_generalinfo
    use m_atomproperties
    use m_optimization
    use m_filenames

    implicit none

    logical exist
    integer, parameter:: nfileLim=10000
    integer, allocatable:: nconf(:)
    integer ifile,IOstatus,nfile,iconf,nottoinclude
    character*80, allocatable:: files(:)
    character*80 line,newstrng

    allocate(files(nfileLim),nconf(nfileLim))

    !Check for filenames starting with 'vasprun' in the directory, and
    !temporarily store them in 'checkfiles.tmp'
    CALL system("ls vasprun*.xml > checkfiles.tmp") 
    inquire(file="checkfiles.tmp",exist=exist)
    if (exist.eqv..false.) then
       print *,'Cannot find vasprun.xml files, stopping.'
       stop
    endif

    !Determine the number of vasprun.xml files to read in
    open(50,file='checkfiles.tmp')
    ifile=1
    do
       read(50,*,IOSTAT=IOstatus) files(ifile)
       if (IOstatus.lt.0) then
          exit
       endif
       if (ifile.gt.nfileLim) then
          print *,'Too many vasp files, please adjust code, stopping.'
          stop
       endif
       ifile=ifile+1
    enddo
    nfile=ifile-1
    close(50)
    CALL system("rm checkfiles.tmp")

    !Determine the number of atomic configurations per file (and record the
    !number that have no configurations: nottoinclude)
    nottoinclude=0
    do ifile=1,nfile
       !print *,'checking file:',trim(files(ifile))
       open(50,file=trim(files(ifile)),IOSTAT=IOstatus)
       iconf=0
       do
         read(50,'(a80)',IOSTAT=IOstatus) line
         if (IOstatus.lt.0) then
            exit
         endif
         newstrng=line(17:22)
         if ( newstrng(1:9).eq.'forces') then
            iconf=iconf+1
         endif
       enddo
       close(50)
       !print *,'file=',ifile,' has ',iconf,' configs'
       if (iconf.eq.0) then
          nottoinclude=nottoinclude+1
       endif
       nconf(ifile)=iconf
    enddo

    !Write the vasprun.xml files to 'fitdbse'
    open(50,file=trim(fitdbse))
    write(50,'(I3,A53)') nfile-nottoinclude,' # Files | Configs to fit | Quantity to fit | Weights'
    do ifile=1,nfile
       !print *,'ifile=',ifile,', nconf=',nconf(ifile)
       if (nconf(ifile).eq.0) then
          !don't include
       elseif (nconf(ifile).lt.10) then
          write(50,'(A30,A3,I1,A12)') files(ifile),' 1-',nconf(ifile),'   Fr      1'
       elseif(nconf(ifile).lt.100) then
          write(50,'(A30,A3,I2,A12)') files(ifile),' 1-',nconf(ifile),'   Fr      1'
       elseif(nconf(ifile).lt.1000) then
          write(50,'(A30,A3,I3,A12)') files(ifile),' 1-',nconf(ifile),'   Fr      1'
       elseif(nconf(ifile).lt.10000) then
          write(50,'(A30,A3,I4,A12)') files(ifile),' 1-',nconf(ifile),'   Fr      1'
       endif
    enddo
    close(50)

    print *
    print *,'Finished writing fitdbse file, stopping.'
    stop

end subroutine createfitdbse
subroutine createMeamParasTemplate

    !--------------------------------------------------------------c
    !
    !     Generates a file 'start_meam_parameters', if requested
    !     by user, which contains a potential in the MEAMfit format
    !     using default values specified in initializemeam and
    !     from an analysis of the separations from the POSCAR 
    !     files.
    !
    !     Called by:     MEAMfit
    !     Calls:         variables_to_p, writemeamp
    !     Returns:       -
    !     Files read:    -
    !     Files written: startparameterfile
    !
    !     Andrew Duff Aug 2014
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_atomproperties
    use m_filenames
    use m_generalinfo
    use m_optimization
    use m_geometry

    implicit none

    integer isp,jsp,i,l,spcCnt,ip
    real(8) smlrad,lrgrad,cutoff

    !These parameters are used to define the 'standard' automaed cut-off generation.
    integer nCutoff     ! Number of cut-off parameters to use for
                                      ! each pair-wise function
    real(8), parameter :: offset = 0.2d0 ! cutoffs need to be made larger than a
                        !given separation in order to affect the enegy at that
                        !separation. This (somewhat arbitrary) constant reflects
                        !that. See below for its use.

    nCutoff=3
    !if (nCutoffOverride.gt.0) then
    !   nCutoff=nCutoffOverride
    !endif

    if (settingsfileexist.eqv..false.) then
       open(unit=1,file=trim(settingsfile),access='append')
       if (lmax.eq.0) then
          !cmin, cmax and meantau not optimized
          freep(1:m3)=0
          freep(1+m3:2*m3)=0
          freep(2*m3+1:2*m3+lm1*m1)=0
       endif
       !electron density parameters:
       freep(2*m3+lm1*m1+1:2*m3+lm1*m1+12*lm1*m2)=0
       spcCnt=0
       do isp=1,m1
          do jsp=1,m1
             if ((isp.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do l=0,lmax
                   do i=1,nCutoff
                      freep(2*m3+lm1*m1+12*(spcCnt*lm1+l)+2*i-1)=2
                      freep(2*m3+lm1*m1+12*(spcCnt*lm1+l)+2*i)=1
                   enddo
                enddo
             endif
             spcCnt=spcCnt+1
          enddo
       enddo
       !one parameter for each embedding function optimized
       freep(2*m3+lm1*m1+12*lm1*m2+1:2*m3+m1 &
           +lm1*m1+12*lm1*m2)=2
       freep(2*m3+m1+lm1*m1+12*lm1*m2+1:2*m3+ &
           2*m1+lm1*m1+12*lm1*m2)=0
       freep(2*m3+2*m1+lm1*m1+12*lm1*m2+1:2*m3+ &
           3*m1+lm1*m1+12*lm1*m2)=0
       freep(2*m3+3*m1+lm1*m1+12*lm1*m2+1:2*m3+ &
           4*m1+lm1*m1+12*lm1*m2)=0
       !pair-potential parameters:
       freep(2*m3+(4+lm1)*m1+12*lm1*m2+1:2*m3+(4+lm1)*m1+32*m2+12*lm1*m2)=0
       spcCnt=0
       do isp=1,m1
          do jsp=1,m1
             if (isp.le.jsp) then
                do i=1,nCutoff
                   freep(2*m3+(4+lm1)*m1+12*lm1*m2+32*spccnt+2*i-1)=2
                   freep(2*m3+(4+lm1)*m1+12*lm1*m2+32*spccnt+2*i)=1
                enddo
             endif
             spcCnt=spcCnt+1
          enddo
       enddo
       if (lmax.eq.0) then
          !rc and rs not optimized
          freep(m4+1:m4+m2)=0
          freep(m4+m2+1:m4+2*m2)=0
       endif
       !enconst randomly generated and optimized
       freep(m4+2*m2+1:m4+2*m2+m1)=2

       write(1,'(I1,A7)') lmax,' ! lmax'
       do i=1,m3
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       do i=1,m3
          write(1,'(I1,A1)',advance='no') freep(i+m3),' '
       enddo
       write(1,*)
       do i=2*m3+1,2*m3+lm1*m1
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       spcCnt=0
       do isp=1,m1
          do jsp=1,m1
             do l=0,lmax
                do i=1,12
                   write(1,'(I1,A1)',advance='no') &
             freep(2*m3+lm1*m1+12*(spcCnt*lm1+l)+i),' '
                enddo
                write(1,*)
             enddo
             spcCnt=spcCnt+1
          enddo
       enddo
       do i=2*m3+lm1*m1+12*lm1*m2+1,2*m3+m1 &
           +lm1*m1+12*lm1*m2
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       do i=2*m3+m1+lm1*m1+12*lm1*m2+1, &
           2*m3+2*m1+lm1*m1+12*lm1*m2
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       do i=2*m3+2*m1+lm1*m1+12*lm1*m2+1, &
           2*m3+3*m1+lm1*m1+12*lm1*m2
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       do i=2*m3+3*m1+lm1*m1+12*lm1*m2+1, &
           2*m3+4*m1+lm1*m1+12*lm1*m2
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       spccnt=0
       do isp=1,m1
          do jsp=1,m1
             do i=1,32
                write(1,'(I1,A1)',advance='no') freep(2*m3+(4+lm1)*m1+ &
                      12*lm1*m2+32*spccnt+i),' '
             enddo
             write(1,*)
             spccnt=spccnt+1
          enddo
       enddo
       !pairpotparameter(16,maxspecies,maxspecies)
       do i=m4+1,m4+m2
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       do i=m4+m2+1,m4+2*m2
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       do i=m4+2*m2+1,m4+2*m2+m1
          write(1,'(I1,A1)',advance='no') freep(i),' '
       enddo
       write(1,*)
       close(1)
       print *
       print *,'Finished writing settings file, stopping.'
       stop
    endif

    print *
    print *,'Setting up meam_parameters file'

    !Set up initial guess at cut-off radii, based on the distribution
    !of separations found from 'findsmallestsepn'
    !---- for the electron densities ----
    if ((thiaccptindepndt.eqv..false.)) then
     print *,' need to generalise code for thiaccptindepndt=false'
     stop
    endif
    do l=0,lmax
       do isp=1,maxspecies
          !First find the smallest separation, note that for the case of
          !electron densities we have to search between separations between
          !species isp and _all_ of the other species.
          smlrad=10000d0
          lrgrad=0d0
          do jsp=isp,maxspecies
             smlrad=min(smlrad,smlsepnperspc_overall(isp,jsp))
             lrgrad=max(lrgrad,lrgsepnperspc_overall(isp,jsp))
          enddo
          do jsp=1,isp
             smlrad=min(smlrad,smlsepnperspc_overall(jsp,isp))
             lrgrad=max(lrgrad,lrgsepnperspc_overall(isp,jsp))
          enddo
          smlrad=smlrad+offset
          lrgrad=min(lrgrad+offset,p_rmax)

          nCutoff=0
          do ip=2*m3+lm1*m1+2+12*l*m2+12*lm1*(isp-1),2*m3+lm1*m1+12+12*l*m2+12*lm1*(isp-1),2
             if (freep(ip).gt.0) then
                nCutoff=nCutoff+1
             endif
          enddo

          i=1
           print *,'for l=',l,', checking:',2*m3+lm1*m1+2+12*l+2*12*lm1*(isp-1),' to ', &
                2*m3+lm1*m1+12+12*l+2*12*lm1*(isp-1)
        ! do ip=2*m3+lm1*m1+2+12*l*m2+12*lm1*(isp-1),2*m3+lm1*m1+12+12*l*m2+12*lm1*(isp-1),2
          do ip=2*m3+lm1*m1+2+12*l+2*12*lm1*(isp-1),2*m3+lm1*m1+12+12*l+2*12*lm1*(isp-1),2
             if (freep(ip).gt.0) then
                cutoff=smlrad + dble(i-1)*(lrgrad-smlrad) / dble(nCutoff-1)
                meamrhodecay(ip-(2*m3+lm1*m1+12*l*m2+12*lm1*(isp-1)),l,1,isp)= cutoff
                if (l.eq.lmax) then
                print *,'Found freep(',ip,')=',freep(ip)
                print *,'setting cutoff=',cutoff,' for meamrhodeacy(', &
      ip-(2*m3+lm1*m1+12*l*m2+12*lm1*(isp-1)),',l=',l,',1,isp=',isp
                endif
                i=i+1
             endif
          enddo
       enddo
    enddo
    !stop
    !--------------------------------------


    !---- for the pair-potentials ----
    spcCnt=0
    do isp=1,maxspecies
        do jsp=1,maxspecies
            if (jsp.ge.isp) then
                smlrad=smlsepnperspc_overall(isp,jsp)+offset
                lrgrad=min(lrgsepnperspc_overall(isp,jsp)+offset,p_rmax)

                nCutoff=0
                do ip=2*m3+(4+lm1)*m1+12*lm1*m2+32*spcCnt+2,2*m3+(4+lm1)*m1+12*lm1*m2+32*(spcCnt+1),2
                   if (freep(ip).gt.0) then
                      nCutoff=nCutoff+1
                   endif
                enddo
                i=1
                do ip=2*m3+(4+lm1)*m1+12*lm1*m2+32*spcCnt+2,2*m3+(4+lm1)*m1+12*lm1*m2+32*(spcCnt+1),2
                   if (freep(ip).gt.0) then
                      cutoff=smlrad + dble(i-1)*(lrgrad-smlrad) / dble(nCutoff-1)
                      pairpotparameter(ip-(2*m3+(4+lm1)*m1+12*lm1*m2+32*spcCnt),isp,jsp)=cutoff
                      i=i+1
                   endif
                enddo
            endif
            spcCnt=spcCnt+1
        enddo
    enddo
    !---------------------------------

    call variables_to_p

    if (writepotfile) then 
       open(2,file=trim(startparameterfile))
       call writemeamp
       close(2)
       print *
       print *,'Written potential file ',trim(startparameterfile),', stopping.'
       stop
    endif

end subroutine createMeamParasTemplate
real(8) function cutoff(r,rmin,rmax)

    !     If rmin<r<rmax, cutoff=1, otherwise cutoff=0
    !
    !     Andrew Duff, 2012

    real(8) r,rmin,rmax

    if ((r.gt.rmin).and.(r.lt.rmax)) then
        cutoff=1d0
    else
        cutoff=0d0
    endif

end function cutoff
subroutine defaultBoundsForRandParas

    !--------------------------------------------------------------c
    !
    !     Initializes default bounds for random generation of
    !     potential parameters.
    !
    !     Called by: MEAMfit
    !     Calls:
    !     Returns:
    !     Files read:
    !     Files written:
    !
    !     Andy Duff, 2014
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_filenames
    use m_atomproperties
    use m_generalinfo
    use m_optimization
    use m_geometry

    implicit none

    integer i,ii,ll,ispc,jspc,spc_cnt,l

    !Variables used for reading in start_meam_parameters file:
    integer nfreep

    allocate( meamtau_minorder(0:lmax,maxspecies), &
        meamtau_maxorder(0:lmax,maxspecies), &
        meamrhodecay_minradius(0:lmax,maxspecies,maxspecies), &
        meamrhodecay_maxradius(0:lmax,maxspecies,maxspecies), &
        meamrhodecay_negvals(0:lmax,maxspecies,maxspecies), &
        meamrhodecay_minorder(1:6,0:lmax,maxspecies,maxspecies), &
        meamrhodecay_maxorder(1:6,0:lmax,maxspecies,maxspecies), &
        meame0_minorder(maxspecies), meame0_maxorder(maxspecies), &
        meamrho0_minorder(maxspecies), meamrho0_maxorder(maxspecies), &
        meamemb3_minorder(maxspecies), meamemb3_maxorder(maxspecies), &
        meamemb4_minorder(maxspecies), meamemb4_maxorder(maxspecies), &
        pairpotparameter_minradius(maxspecies,maxspecies), &
        pairpotparameter_maxradius(maxspecies,maxspecies), &
        pairpotparameter_negvals(maxspecies,maxspecies), &
        pairpotparameter_minorder(1:8,maxspecies,maxspecies), &
        pairpotparameter_maxorder(1:8,maxspecies,maxspecies), &
        minordervaluepairpot(1:16), &
        maxordervaluepairpot(1:16), &
        nfreeppairpot(maxspecies,maxspecies), &
        enconst_minorder(maxspecies), enconst_maxorder(maxspecies) )

    !Determine max and min values for cutoffs based on separation distribution
    !from fitting database
    meamrhodecay_maxradius=0d0
    pairpotparameter_maxradius=0d0
    do ispc=1,maxspecies
       !First, max values
       meamrhodecay_maxradius(0,1,ispc)=0d0
       do jspc=ispc,maxspecies
          meamrhodecay_maxradius(0,1,ispc)=max(meamrhodecay_maxradius(0,1,ispc),lrgsepnperspc_overall(ispc,jspc))
       enddo
       do jspc=1,ispc
          meamrhodecay_maxradius(0,1,ispc)=max(meamrhodecay_maxradius(0,1,ispc),lrgsepnperspc_overall(ispc,jspc))
       enddo
       !Now, min values
       meamrhodecay_minradius(0,1,ispc)=10000d0
       do jspc=ispc,maxspecies
          meamrhodecay_minradius(0,1,ispc)=min(meamrhodecay_minradius(0,1,ispc),smlsepnperspc_overall(ispc,jspc))
       enddo
       do jspc=1,ispc
          meamrhodecay_minradius(0,1,ispc)=min(meamrhodecay_minradius(0,1,ispc),smlsepnperspc_overall(jspc,ispc))
       enddo
       !These min values can't be below 'cutoffMin' (specified either by user or
       !given a default value - neccessary to stop cutoffs enchroaching on
       !Biersack-Ziegler low separation limit.
       if (meamrhodecay_minradius(0,1,ispc).lt.cutoffMin) then
          meamrhodecay_minradius(0,1,ispc)=cutoffMin
       endif

    enddo
    do l=1,lmax
       do ispc=1,maxspecies
          meamrhodecay_minradius(l,1,ispc)=meamrhodecay_minradius(0,1,ispc)
          meamrhodecay_maxradius(l,1,ispc)=meamrhodecay_maxradius(0,1,ispc)
       enddo
    enddo
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                pairpotparameter_minradius(ispc,jspc)=max(smlsepnperspc_overall(ispc,jspc), &
                                                          cutoffMin)
                pairpotparameter_maxradius(ispc,jspc)=lrgsepnperspc_overall(ispc,jspc)
            endif
        enddo
    enddo

    !meamtau
    if (lmax.gt.0) then
       do ll=1,lmax
           do ispc=1,maxspecies
               if (freep(2*m3+ispc+ll*m1).ge.2) then
                   meamtau_minorder(ll,ispc)=meamtau_minorder_default
                   meamtau_maxorder(ll,ispc)=meamtau_maxorder_default
               endif
           enddo
       enddo
    endif

    !electron density functions
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    !Check if we need to generate random paramaters for
                    !thi(ispc,jspc,ll)
                    nfreep=0
                    do i=1,12
                        if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i).ge.2) then
                            nfreep=nfreep+1
                        endif
                    enddo
                    if (nfreep.gt.0) then
                        if (typethi(ispc,jspc).eq.1) then !Cubic form
                            meamrhodecay_negvals(ll,ispc,jspc)=meamrhodecay_negvals_default

                            !For each coefficient which needs to be generated,
                            !read in
                            !min and max orders (10^..)
                            ii=1
                            do i=1,6
                                if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+2*i-1).ge.2) &
                                    then
                                    meamrhodecay_minorder(ii,ll,ispc,jspc)=meamrhodecay_minorder_default
                                    meamrhodecay_maxorder(ii,ll,ispc,jspc)=meamrhodecay_maxorder_default
                                    ii=ii+1
                                endif
                            enddo
                        endif
                    endif
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !embedding functions
    embfuncRandType=2 !Logarithmic by default
    do ispc=1,maxspecies
        do i=1,4
            if (freep(2*m3+lm1*m1+12*lm1*m2+ispc+(i-1)*m1).ge.2) then
                if (i.eq.1) then
                    meame0_minorder(ispc)=meame0_minorder_default
                    meame0_maxorder(ispc)=meame0_maxorder_default
                elseif (i.eq.2) then
                    meamrho0_minorder(ispc)=meamrho0_minorder_default
                    meamrho0_maxorder(ispc)=meamrho0_maxorder_default
                elseif (i.eq.3) then
                    meamemb3_minorder(ispc)=meamemb3_minorder_default
                    meamemb3_maxorder(ispc)=meamemb3_maxorder_default
                elseif (i.eq.4) then
                    meamemb4_minorder(ispc)=meamemb4_minorder_default
                    meamemb4_maxorder(ispc)=meamemb4_maxorder_default
                endif
            endif
        enddo
    enddo

    !pair-potentials
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                !Check if we need to generate random paramaters for
                !pairpot(isp,jspc)
                nfreeppairpot(ispc,jspc)=0
                do i=1,32
                    if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt) &
                        .ge.2) then
                        nfreeppairpot(ispc,jspc)=nfreeppairpot(ispc,jspc)+1
                    endif
                enddo
                if (nfreeppairpot(ispc,jspc).gt.0) then
                    if (typepairpot(ispc,jspc).eq.2) then
                        pairpotparameter_negvals(ispc,jspc)=pairpotparameter_negvals_default
                        if ((pairpotparameter_negvals(ispc,jspc).lt.1).or. &
                            (pairpotparameter_negvals(ispc,jspc).gt.3)) then
                            print *,'      ERROR: you must specify 1, 2 or 3,'
                            print *,'      stopping.'
                            stop
                        endif
                        !Setup default minimum and maximum orders for generating
                        !pair-potential coefficients
                        ii=1
                        do i=1,16
                            if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+2*i-1+32*spc_cnt) &
                                .ge.2) then
                                pairpotparameter_minorder(ii,ispc,jspc)=pairpotparameter_minorder_default
                                pairpotparameter_maxorder(ii,ispc,jspc)=pairpotparameter_maxorder_default
                                ii=ii+1
                            endif
                        enddo
                    else
                        print *,'typepairpot=',typepairpot(ispc,jspc), &
                            'not supported,stopping'
                        stop
                    endif
                endif
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !enconst
    do ispc=1,maxspecies
        if (freep(m4+2*m2+ispc).ge.2) then
            enconst_minorder(ispc)=enconst_minorder_default
            enconst_maxorder(ispc)=enconst_maxorder_default
        endif
    enddo

end subroutine defaultBoundsForRandParas
subroutine displayBoundsForRandParas

    !--------------------------------------------------------------c
    !
    !     Write to standard output the limits which are to be used 
    !     to initialize the random starting values of the MEAM 
    !     parameters.
    !
    !     Called by: MEAMfit
    !     Calls: -
    !     Returns: -
    !     Files read: -
    !     Files written: -
    !
    !     Andrew Duff, 2014
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_filenames
    use m_atomproperties
    use m_generalinfo
    use m_optimization

    implicit none

    logical meamtauLimitsSame,meamtauRandParasSame
    integer i,j,ii,ll,ispc,jspc,spc_cnt,minorderStr,maxorderStr
    integer nfreep

    !cmin and cmax
    j=1
    ispc=1
    jspc=1
    do i=1,m3
        if (freep(i).ge.2) then
            write(*,*) 'cmin(',j,',',ispc,',',jspc,') ', &
                'to be randomly seeded'
        endif
        if (j.lt.m1) then
            j=j+1
        else
            if (ispc.lt.m1) then
                j=1
                ispc=ispc+1
            else
                j=1
                ispc=1
                jspc=jspc+1
            endif
        endif
    enddo
    j=1
    ispc=1
    jspc=1
    do i=m3+1,2*m3
        if (freep(i).ge.2) then
            write(*,*) 'cmax(',j,',',ispc,',',jspc,') ', &
                'to be randomly seeded'
        endif   
        if (j.lt.m1) then
            j=j+1   
        else        
            if (ispc.lt.m1) then
                j=1 
                ispc=ispc+1
            else
                j=1
                ispc=1
                jspc=jspc+1
            endif
        endif
    enddo

    !meamtau

    !First check if all meamtau paras are to be randonly generated and that
    !their limits for random para gen are the same.
    if (lmax.gt.0) then

       meamtauRandParasSame=.true.
       meamtauLimitsSame=.true.
       minorderStr=meamtau_minorder(1,1)
       maxorderStr=meamtau_maxorder(1,1)
       do ll=1,lmax
          do ispc=1,maxspecies
             if (freep(2*m3+ispc+ll*m1).lt.2) then
                meamtauRandParasSame=.false.
             endif
             if ((meamtau_minorder(ll,ispc).ne.minorderStr).or. &
                 (meamtau_maxorder(ll,ispc).ne.maxorderStr)) then
                meamtauLimitsSame=.false.
             endif
          enddo
       enddo

       if (meamtauRandParasSame.and.meamtauLimitsSame) then
          write(*,'(A19,I1,A5,I1,A1)') '    meamtau(ispc=1-',maxspecies,',l=1-',lmax,')'
          write(*,'(A25,I2,A6,I2)') '       Can vary from: 10^', &
             meamtau_minorder(1,1),' -10^',meamtau_maxorder(1,1)
       else
          do ll=1,lmax
              do ispc=1,maxspecies
                  if (freep(2*m3+ispc+ll*m1).ge.2) then
                      print *,'    meamtau(ispc=',ispc,',l=',ll,')'
                      write(*,'(A25,I2,A6,I2)') '       Can vary from: 10^', &
                          meamtau_minorder(ll,ispc),' -10^',meamtau_maxorder(ll,ispc)
                  endif
              enddo
          enddo
       endif

    endif

    !electron density functions
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    !Check if we need to generate random paramaters for
                    !thi(ispc,jspc,ll)
                    nfreep=0
                    do i=1,12
                        if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i).ge.2) then
                            nfreep=nfreep+1
                        endif
                    enddo
                    if (nfreep.gt.0) then
                        write(*,'(A8,I1,A1,I1,A3,I1,A2)') &
                            '    thi(',ispc,',',jspc,',l=',ll,'):'
                        if (typethi(ispc,jspc).eq.1) then !Cubic form
                            !Read in parameters from 'dataforMEAMparagen'
                            !necessary
                            !to generate these parameters
                            write(*,*) '      Function: Sum over cubic terms'
                            write(*,'(A28,I2)') '       Number of free paras:',nfreep
                            write(*,'(A36,F5.3,A1,F5.3)') &
                                '       Allowed cutoff radius range: ', &
                                meamrhodecay_minradius(ll,ispc,jspc),'-', &
                                meamrhodecay_maxradius(ll,ispc,jspc)
                            if (meamrhodecay_negvals(ll,ispc,jspc).eq.1) then
                                write(*,*) '      No negative values allowed.'
                            elseif (meamrhodecay_negvals(ll,ispc,jspc).eq.2) then
                                write(*,*) '      Only negative values allowed.'
                            elseif (meamrhodecay_negvals(ll,ispc,jspc).eq.3) then
                                write(*,*) '      Both positive and negative ', &
                                    'values allowed.'
                            endif

                            !For each coefficient which needs to be generated,
                            !read in
                            !min and max orders (10^..)
                            ii=1
                            do i=1,6
                                if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+2*i-1).ge.2) &
                                    then
                                    write(*,'(A19,I2,A19,I2,A6,I2)') &
                                        '       Coefficient ',ii,' can vary from: 10^', &
                                        meamrhodecay_minorder(ii,ll,ispc,jspc),' - 10^', &
                                        meamrhodecay_maxorder(ii,ll,ispc,jspc)
                                    ii=ii+1
                                endif
                            enddo
                        endif
                    endif
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !embedding functions
    do ispc=1,maxspecies
        do i=1,4
            if (freep(2*m3+lm1*m1+12*lm1*m2+ispc+(i-1)*m1).ge.2) then
                if (embfuncRandType.eq.1) then
                 if (i.eq.1) then
                     write(*,'(A10,I1,A2)') '    E_emb(',ispc,'):'
                     write(*,'(A29,F10.5,A3,F10.5)') &
                         '       meame0 can vary from: ', &
                         meame0_minorder(ispc),' - ', &
                         meame0_maxorder(ispc)
                 elseif (i.eq.2) then
                     write(*,'(A31,F10.5,A3,F10.5)') &
                         '       meamrho0 can vary from: ', &
                         meamrho0_minorder(ispc),' - ', &
                         meamrho0_maxorder(ispc)
                 elseif (i.eq.3) then
                     write(*,'(A31,F10.5,A3,F10.5)') &
                         '       meamemb3 can vary from: ', &
                         meamemb3_minorder(ispc),' - ', &
                         meamemb3_maxorder(ispc)
                 elseif (i.eq.4) then
                     write(*,'(A31,F10.5,A3,F10.5)') &
                         '       meamemb4 can vary from: ', &
                         meamemb4_minorder(ispc),' - ', &
                         meamemb4_maxorder(ispc)
                 endif
                else
                 if (i.eq.1) then
                     write(*,'(A10,I1,A2)') '    E_emb(',ispc,'):'
                     write(*,'(A32,F10.5,A6,F10.5)') &
                         '       meame0 can vary from: 10^', &
                         meame0_minorder(ispc),' - 10^', &
                         meame0_maxorder(ispc)
                 elseif (i.eq.2) then
                     write(*,'(A34,F10.5,A6,F10.5)') &
                         '       meamrho0 can vary from: 10^', &
                         meamrho0_minorder(ispc),' - 10^', &
                         meamrho0_maxorder(ispc)
                 elseif (i.eq.3) then
                     write(*,'(A34,F10.5,A6,F10.5)') &
                         '       meamemb3 can vary from: 10^', &
                         meamemb3_minorder(ispc),' - 10^', &
                         meamemb3_maxorder(ispc)
                 elseif (i.eq.4) then
                     write(*,'(A34,F10.5,A6,F10.5)') &
                         '       meamemb4 can vary from: 10^', &
                         meamemb4_minorder(ispc),' - 10^', &
                         meamemb4_maxorder(ispc)
                 endif
                endif
            endif
        enddo
    enddo

    !pair-potentials
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                !Check if we need to generate random paramaters for
                !pairpot(isp,jspc)
                nfreeppairpot(ispc,jspc)=0
                do i=1,32
                    if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt) &
                        .ge.2) then
                        nfreeppairpot(ispc,jspc)=nfreeppairpot(ispc,jspc)+1
                    endif
                enddo
                if (nfreeppairpot(ispc,jspc).gt.0) then
                    write(*,'(A6,I1,A1,I1,A2)') &
                        '    V(',ispc,',',jspc,'):'
                    if (typepairpot(ispc,jspc).eq.2) then
                        !Read in parameters from 'dataforMEAMparagen' necessary
                        !to
                        !generate these parameters
                        write(*,*) '      Function: Sum over cubic terms'
                        write(*,'(A28,I2)') '       Number of free paras:', &
                            nfreeppairpot(ispc,jspc)
                        write(*,'(A36,F5.3,A1,F5.3)') &
                            '       Allowed cutoff radius range: ', &
                            pairpotparameter_minradius(ispc,jspc),'-', &
                            pairpotparameter_maxradius(ispc,jspc)
                        if (pairpotparameter_negvals(ispc,jspc).eq.1) then
                            write(*,*) '      No negative values allowed.'
                        elseif (pairpotparameter_negvals(ispc,jspc).eq.2) then
                            write(*,*) '      Only negative values allowed.'
                        elseif (pairpotparameter_negvals(ispc,jspc).eq.3) then
                            write(*,*) '      Both positive and negative ', &
                                'values allowed.'
                        else
                            print *,'      ERROR: you must specify 1, 2 or 3,'
                            print *,'      stopping.'
                            stop
                        endif
                        !For each coefficient which needs to be generated, read
                        !in min and max orders (10^..)
                        ii=1
                        do i=1,16
                            if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+2*i-1+32*spc_cnt) &
                                .ge.2) then
                                write(*,'(A19,I2,A19,I2,A6,I2)') &
                                    '       Coefficient ',ii,' can vary from: 10^', &
                                    pairpotparameter_minorder(ii,ispc,jspc),' - 10^', &
                                    pairpotparameter_maxorder(ii,ispc,jspc)
                                ii=ii+1
                            endif
                        enddo
                    else
                        print *,'typepairpot=',typepairpot(ispc,jspc), &
                            'not supported,stopping'
                        stop
                    endif
                endif
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !enconst
    do ispc=1,maxspecies
        if (freep(m4+2*m2+ispc).ge.2) then
            write(*,'(A12,I1,A20,F6.3,A6,F6.3)') &
                '    enconst(',ispc,') can vary from: 10^', &
                enconst_minorder(ispc),' - 10^',enconst_maxorder(ispc)
        endif
    enddo

end subroutine displayBoundsForRandParas
subroutine displaytime(unit,hrs)

    !-------------------------------------------------------------------c
    !
    !     Displays total time for code to finish.
    !
    !     Called by:     program MEAMfit
    !     Calls:         -
    !     Arguments:     -
    !     Returns:       -
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff 2014
    !
    !-------------------------------------------------------------------c

    use m_generalinfo

    implicit none

    integer unit,tottime,secs,mins,hrs

    call getTime(tottime)
    tottime=tottime-tstart
    ! 
    ! call CPU_TIME(tottime)

    !Display different units for the total time depending on whether only
    !seconds, minutes or hours have already passed.
    secs=0
    mins=0
    hrs=0
    if (tottime.lt.60) then
       write(unit,'(A18,I2,A2)') ' Total time taken ',INT(tottime),'s.'
    elseif (tottime.lt.3600) then
       mins=INT(tottime/60)
       secs=MOD(tottime,60)
       write(unit,'(A18,I2,A2,I2,A2)') ' Total time taken ',INT(mins),'m ', &
                 INT(secs),'s.'
    elseif (tottime.ge.3600) then
       hrs=INT(tottime/3600)
       mins=MOD(tottime,3600)/60d0
       secs=MOD(mins,1)*60d0
       write(unit,'(A18,I5,A2,I2,A2,I2,A2)') ' Total time taken ',INT(hrs),'h ', &
                 INT(mins),'m ',INT(secs),'s.'
    endif

end subroutine displaytime
subroutine electrondensity(iatom,cart)

    !----------------------------------------------------------------c
    !
    !     Fills the module m_electrondensity with the rho_i_l
    !     parameters and also computes the meam t s. The code is
    !     hardwired for Lmax up to and including 4.
    !
    !     Called by:     meamenergy
    !     Calls:         distance, distance2
    !     Returns:       rhol,meam_t
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter, Feb 9 2006
    !
    !----------------------------------------------------------------c

    use m_generalinfo
    use m_atomproperties
    use m_geometry
    use m_meamparameters
    use m_screening       !S_ij
    use m_electrondensity !f_i_L_ij,rho_i_l

    implicit none

    integer i,j,jj,nni,iaux(1),isp,jsp,iatom,cart
    real(8) aux1(3),aux2(3,3),aux3(3,3,3),aux4(3,3,3,3), &
        aux2a,aux3a(3),aux,tmp1,tmp2

    if(allocated(rhol)) deallocate(rhol)
    iaux=maxval(gn_neighbors(:,istr))
    allocate(rhol(0:lmax,gn_inequivalentsites(istr)), &
        meam_t(1:lmax,gn_inequivalentsites(istr)))
    do i=1,gn_inequivalentsites(istr)

        isp=gspecies(i,istr)
        nni=gn_neighbors(i,istr)
        rhol(0:lmax,i)=0d0
        meam_t(1:lmax,i)=0d0   !angmom=0 never used
        aux1=0d0
        aux2=0d0
        aux3=0d0
        aux4=0d0
        aux2a=0d0
        aux3a=0d0

        !Only use stored separations _if_ considering undisplaced structures
        !_OR_ fastForce set to .true. (separations stored also for displaced
        !structures)
        if ( ((iatom.eq.0).and.(cart.eq.0)) .or. fastForce ) then

           do jj=1,nni
               j=gneighborlist(jj,i,istr)
               jsp=gspecies(j,istr)
               rij=diststr(jj,i,iatom,cart,istr)

               !          if (rij.ne.diststr(jj,i,iatom,cart,istr)) then
               !             print *,'rij does not match the stored value'
               !             print *,'jj=',jj,', i=',i,', istr=',istr,
               !     + ', iatom=',iatom,',cart=',cart
               !             write(*,'(A9,F20.15,A9,F20.15)') 'rij calc=',rij,
               !     +         ', stored=',diststr(jj,i,iatom,cart,istr)
               !             stop
               !          endif

               rij2=dist2str(jj,i,iatom,cart,istr)

               !          if (rij2.ne.dist2str(jj,i,iatom,cart,istr)) then
               !             print *,'rij2 does not match the stored value'
               !             print *,'jj=',jj,', i=',i,', istr=',istr,
               !     + ', iatom=',iatom,',cart=',cart
               !             write(*,'(A9,F20.15,A9,F20.15)') 'rij calc=',rij2,
               !     +         ', stored=',dist2str(jj,i,iatom,cart,istr)
               !             stop
               !          endif

               rij3=dist3str(jj,i,iatom,cart,istr)

               rhol(0,i)=rhol(0,i)+screening(jj,i)*fjlij(0,jj,i)

               if (lmax.ge.1) then

                   if ((gxyz(1,j,istr)-gxyz(1,i,istr)).ne. &
                       dxstr(jj,i,iatom,cart,istr)) then
                       print *,'(gxyz(1,j,istr)-gxyz(1,i,istr).ne.', &
                           'dxstr(jj,i,iatom,cart,istr), stopping.'
                       print *,'j=',j,', i=',i,', jj=',jj,', istr=',istr, &
                           ', iatom=',iatom,', cart=',cart
                       stop
                   endif
                   aux1(1)=aux1(1)+screening(jj,i)*fjlij(1,jj,i)/rij* &
                       dxstr(jj,i,iatom,cart,istr)
                   aux1(2)=aux1(2)+screening(jj,i)*fjlij(1,jj,i)/rij* &
                       dystr(jj,i,iatom,cart,istr)
                   aux1(3)=aux1(3)+screening(jj,i)*fjlij(1,jj,i)/rij* &
                       dzstr(jj,i,iatom,cart,istr)

                   meam_t(1,i)=meam_t(1,i)+meamtau(1,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
                   if (rij.eq.0) then
                       print *,'found an rij equal to zero'
                       print *,i,gneighborlist(jj,i,1)
                       stop
                   endif
               endif
               if (lmax.ge.2) then

                   aux2(1,1)=aux2(1,1)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
                       dxstr(jj,i,iatom,cart,istr)*dxstr(jj,i,iatom,cart,istr)
                   aux2(1,2)=aux2(1,2)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
                       dxstr(jj,i,iatom,cart,istr)*dystr(jj,i,iatom,cart,istr)
                   aux2(1,3)=aux2(1,3)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
                       dxstr(jj,i,iatom,cart,istr)*dzstr(jj,i,iatom,cart,istr)
                   aux2(2,2)=aux2(2,2)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
                       dystr(jj,i,iatom,cart,istr)*dystr(jj,i,iatom,cart,istr)
                   aux2(2,3)=aux2(2,3)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
                       dystr(jj,i,iatom,cart,istr)*dzstr(jj,i,iatom,cart,istr)
                   aux2(3,3)=aux2(3,3)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
                       dzstr(jj,i,iatom,cart,istr)*dzstr(jj,i,iatom,cart,istr)

                   aux2a=aux2a+screening(jj,i)*fjlij(2,jj,i)
                   meam_t(2,i)=meam_t(2,i)+meamtau(2,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
               endif
               if (lmax.ge.3) then

                   aux3(1,1,1)=aux3(1,1,1)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dxstr(jj,i,iatom,cart,istr)*dxstr(jj,i,iatom,cart,istr)* &
                       dxstr(jj,i,iatom,cart,istr)
                   aux3(1,1,2)=aux3(1,1,2)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dxstr(jj,i,iatom,cart,istr)*dxstr(jj,i,iatom,cart,istr)* &
                       dystr(jj,i,iatom,cart,istr)
                   aux3(1,1,3)=aux3(1,1,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dxstr(jj,i,iatom,cart,istr)*dxstr(jj,i,iatom,cart,istr)* &
                       dzstr(jj,i,iatom,cart,istr)
                   aux3(1,2,2)=aux3(1,2,2)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dxstr(jj,i,iatom,cart,istr)*dystr(jj,i,iatom,cart,istr)* &
                       dystr(jj,i,iatom,cart,istr)
                   aux3(1,2,3)=aux3(1,2,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dxstr(jj,i,iatom,cart,istr)*dystr(jj,i,iatom,cart,istr)* &
                       dzstr(jj,i,iatom,cart,istr)
                   aux3(1,3,3)=aux3(1,3,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dxstr(jj,i,iatom,cart,istr)*dzstr(jj,i,iatom,cart,istr)* &
                       dzstr(jj,i,iatom,cart,istr)
                   aux3(2,2,2)=aux3(2,2,2)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dystr(jj,i,iatom,cart,istr)*dystr(jj,i,iatom,cart,istr)* &
                       dystr(jj,i,iatom,cart,istr)
                   aux3(2,2,3)=aux3(2,2,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dystr(jj,i,iatom,cart,istr)*dystr(jj,i,iatom,cart,istr)* &
                       dzstr(jj,i,iatom,cart,istr)
                   aux3(2,3,3)=aux3(2,3,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dystr(jj,i,iatom,cart,istr)*dzstr(jj,i,iatom,cart,istr)* &
                       dzstr(jj,i,iatom,cart,istr)
                   aux3(3,3,3)=aux3(3,3,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
                       dzstr(jj,i,iatom,cart,istr)*dzstr(jj,i,iatom,cart,istr)* &
                       dzstr(jj,i,iatom,cart,istr)
                   aux3a(1)=aux3a(1)+screening(jj,i)*fjlij(3,jj,i)/rij* &
                       dxstr(jj,i,iatom,cart,istr)
                   aux3a(2)=aux3a(2)+screening(jj,i)*fjlij(3,jj,i)/rij* &
                       dystr(jj,i,iatom,cart,istr)

                   meam_t(3,i)=meam_t(3,i)+meamtau(3,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
               endif
               if (lmax.ge.4) then
                   print *,'not yet!'
                   meam_t(4,i)=meam_t(4,i)+meamtau(4,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
                   stop
               endif
           enddo  !sum over neighbors complete

        else

           do jj=1,nni

               j=gneighborlist(jj,i,istr)
               jsp=gspecies(j,istr)

               rij=distance(i,j)
               rij2=distance2(i,j)
               rij3=rij2*rij

               rhol(0,i)=rhol(0,i)+screening(jj,i)*fjlij(0,jj,i)

               if (lmax.ge.1) then

                   aux1(1)=aux1(1)+screening(jj,i)*fjlij(1,jj,i)/rij* &
                       (gxyz(1,j,istr)-gxyz(1,i,istr))
                   aux1(2)=aux1(2)+screening(jj,i)*fjlij(1,jj,i)/rij* &
                       (gxyz(2,j,istr)-gxyz(2,i,istr))
                   aux1(3)=aux1(3)+screening(jj,i)*fjlij(1,jj,i)/rij* &
                       (gxyz(3,j,istr)-gxyz(3,i,istr))

                   meam_t(1,i)=meam_t(1,i)+meamtau(1,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
                   if (rij.eq.0) then
                       print *,'found an rij equal to zero'
                       print *,i,gneighborlist(jj,i,1)
                       stop
                   endif
               endif
               if (lmax.ge.2) then

                   aux2(1,1)=aux2(1,1)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(1,j,istr)-gxyz(1,i,istr))
                   aux2(1,2)=aux2(1,2)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(2,j,istr)-gxyz(2,i,istr))
                   aux2(1,3)=aux2(1,3)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(3,j,istr)-gxyz(3,i,istr))
                   aux2(2,2)=aux2(2,2)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
              (gxyz(2,j,istr)-gxyz(2,i,istr))*(gxyz(2,j,istr)-gxyz(2,i,istr))
                   aux2(2,3)=aux2(2,3)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
              (gxyz(2,j,istr)-gxyz(2,i,istr))*(gxyz(3,j,istr)-gxyz(3,i,istr))
                   aux2(3,3)=aux2(3,3)+screening(jj,i)*fjlij(2,jj,i)/rij2* &
              (gxyz(3,j,istr)-gxyz(3,i,istr))*(gxyz(3,j,istr)-gxyz(3,i,istr))

                   aux2a=aux2a+screening(jj,i)*fjlij(2,jj,i)
                   meam_t(2,i)=meam_t(2,i)+meamtau(2,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
               endif
               if (lmax.ge.3) then

                  aux3(1,1,1)=aux3(1,1,1)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(1,j,istr)-gxyz(1,i,istr)) &
             *(gxyz(1,j,istr)-gxyz(1,i,istr))
                  aux3(1,1,2)=aux3(1,1,2)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(1,j,istr)-gxyz(1,i,istr)) &
             *(gxyz(2,j,istr)-gxyz(2,i,istr))
                  aux3(1,1,3)=aux3(1,1,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(1,j,istr)-gxyz(1,i,istr)) &
             *(gxyz(3,j,istr)-gxyz(3,i,istr))
                  aux3(1,2,2)=aux3(1,2,2)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(2,j,istr)-gxyz(2,i,istr)) &
             *(gxyz(2,j,istr)-gxyz(2,i,istr))
                  aux3(1,2,3)=aux3(1,2,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(2,j,istr)-gxyz(2,i,istr)) &
             *(gxyz(3,j,istr)-gxyz(3,i,istr))
                  aux3(1,3,3)=aux3(1,3,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))*(gxyz(3,j,istr)-gxyz(3,i,istr)) &
             *(gxyz(3,j,istr)-gxyz(3,i,istr))
                  aux3(2,2,2)=aux3(2,2,2)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(2,j,istr)-gxyz(2,i,istr))*(gxyz(2,j,istr)-gxyz(2,i,istr)) &
             *(gxyz(2,j,istr)-gxyz(2,i,istr))
                  aux3(2,2,3)=aux3(2,2,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(2,j,istr)-gxyz(2,i,istr))*(gxyz(2,j,istr)-gxyz(2,i,istr)) &
             *(gxyz(3,j,istr)-gxyz(3,i,istr))
                  aux3(2,3,3)=aux3(2,3,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(2,j,istr)-gxyz(2,i,istr))*(gxyz(3,j,istr)-gxyz(3,i,istr)) &
             *(gxyz(3,j,istr)-gxyz(3,i,istr))
                  aux3(3,3,3)=aux3(3,3,3)+screening(jj,i)*fjlij(3,jj,i)/rij3* &
              (gxyz(3,j,istr)-gxyz(3,i,istr))*(gxyz(3,j,istr)-gxyz(3,i,istr)) &
             *(gxyz(3,j,istr)-gxyz(3,i,istr))
                  aux3a(1)=aux3a(1)+screening(jj,i)*fjlij(3,jj,i)/rij* &
              (gxyz(1,j,istr)-gxyz(1,i,istr))
                  aux3a(2)=aux3a(2)+screening(jj,i)*fjlij(3,jj,i)/rij* &
              (gxyz(2,j,istr)-gxyz(2,i,istr))

                   meam_t(3,i)=meam_t(3,i)+meamtau(3,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
               endif
               if (lmax.ge.4) then
                   print *,'not yet!'
                   meam_t(4,i)=meam_t(4,i)+meamtau(4,jsp)* &
                       screening(jj,i)*fjlij(0,jj,i)
                   stop
               endif
           enddo  !sum over neighbors complete

        endif

        if(lmax.ge.1)then
            if (rhol(0,i).ne.0d0) then
                meam_t(1,i)=meam_t(1,i)/rhol(0,i)
            endif
            rhol(1,i)=sqrt(aux1(1)**2+aux1(2)**2+aux1(3)**2)
        endif
        if(lmax.ge.2)then
            if (rhol(0,i).ne.0d0) then
                meam_t(2,i)=meam_t(2,i)/rhol(0,i)
            endif
            tmp1=aux2(1,1)**2+2d0*aux2(1,2)**2+2d0*aux2(1,3)**2+ &
                aux2(2,2)**2+2d0*aux2(2,3)**2+aux2(3,3)**2
            tmp2=aux2a*aux2a/3d0
            if (abs((tmp1-tmp2)/tmp1).lt.10d-14) then !first check that
                aux=0d0           !tmp1 and tmp2 aren't equal (up to the
            else                 !15d.p. prec of dble prec). If they are,
                aux=tmp1-tmp2     !then have to set tmp1-tmp2 equal to 0
            endif                !by hand (otherwise it will be of O
            !10d-16, which will incorr. give sqr rt=
            !10d-8)
            rhol(2,i)=sqrt(abs(aux)) !I use BASKES formula 8c PRB 46,
            !2727 (1992)

            !           if ((i.eq.1).and.(istr.eq.52)) then
            !               print *,'rhol(2,1)^2=',rhol(2,i)**2
            !           endif

        endif
        if(lmax.ge.3)then
            if (rhol(0,i).ne.0d0) then
                meam_t(3,i)=meam_t(3,i)/rhol(0,i)
            endif
            aux=aux3(1,1,1)**2+3d0*aux3(1,1,2)**2+3d0*aux3(1,1,3)**2+ &
                3d0*aux3(1,2,2)**2+6d0*aux3(1,2,3)**2+3d0*aux3(1,3,3)**2+ &
                aux3(2,2,2)**2+3d0*aux3(2,2,3)**2+3d0*aux3(2,3,3)**2+ &
                aux3(3,3,3)**2
            rhol(3,i)=sqrt(abs(aux))    !I use BASKES formula 8d PRB 46,
            !2727 (1992)
        endif
        if(lmax.ge.4)then
            if (rhol(0,i).ne.0d0) then
                meam_t(4,i)=meam_t(4,i)/rhol(0,i)
            endif
            print *,'not yet!'
            stop
        endif
    enddo

end subroutine electrondensity
subroutine embeddingfunction

    !-------------------------------------------------------------c
    !
    !     Calculates the embedding-function, F(rho_i) (=meam_f(i)
    !     in the code), for each of the inequivalent sites for
    !     of the 'istr'th structure.
    !     (with Camelion modification).
    !
    !     Called by:     meamenergy
    !     Calls:         -
    !     Returns:       meam_f
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter and Andrew Duff, 2006-2015
    !
    !-------------------------------------------------------------c

    use m_generalinfo
    use m_atomproperties
    use m_geometry
    use m_meamparameters
    use m_electrondensity

    implicit none

    integer i,isp

    if (embfunctype.eq.1) then

        do i=1,gn_inequivalentsites(istr)
            isp=gspecies(i,istr)
            meam_f(i)= -meame0(isp)*sqrt(rho_i(i)) &
                + meamrho0(isp)*(rho_i(i)**2) &
                + meamemb3(isp)*(rho_i(i)**3) &
                + meamemb4(isp)*(rho_i(i)**4)
        enddo
        
    elseif (embfunctype.eq.2) then

        !The Hepburn embedding function

        !Included during early stage of code development - can re-instate later
        !if there is demand but for now not sufficiently tested to include.
        print *,'ERROR: embedding function type unsupported, stopping.'
        stop

        do i=1,gn_inequivalentsites(istr)

            isp=gspecies(i,istr)
            if (isp.eq.1) then

                meam_f(i)=-meame0(isp)*sqrt(rho_i(i)) &
                    + meamrho0(isp)*(rho_i(i)**2) &
                    + meamemb3(isp)*(rho_i(i)**4)
            elseif (isp.eq.2) then

                if (rho_i(i).le.0d0) then
                    meam_f(i)=0d0
                elseif ((rho_i(i).gt.0d0).and. &
                        (rho_i(i).le.0.001d0)) then
                    meam_f(i)= -meame0(2)* ( &
                        3d0*((rho_i(i)/0.001d0)**2)- &
                        2d0*((rho_i(i)/0.001d0)**3) )
                elseif ((rho_i(i).gt.0.001d0).and. &
                        (rho_i(i).le.0.5)) then
                    meam_f(i)= -meame0(2)
                elseif ((rho_i(i).gt.0.5d0).and. &
                        (rho_i(i).le.1d0)) then
                    meam_f(i)= 5d0*(-meamrho0(2)) - meame0(2) &
                        - 12d0*(-meamrho0(2))*(rho_i(i)/0.5d0) &
                        + 9d0*(-meamrho0(2))*( (rho_i(i)/0.5d0)**2 ) &
                        - 2d0*(-meamrho0(2))*( (rho_i(i)/0.5d0)**3 )
                elseif (rho_i(i).gt.1d0) then
                    meam_f(i)=-meame0(2)-meamrho0(2)
                endif
            else
                print *,'isp.ne.1 or 2 in emb en routine'
                stop
            endif

        enddo

    elseif( embfunctype.eq.3) then

        !Alternative embedding function formalism from Baskes

        !Included during early stage of code development - can re-instate later
        !if there is demand but for now not sufficiently tested to include.
        print *,'ERROR: embedding function type unsupported, stopping.'
        stop

        do i=1,gn_inequivalentsites(istr)

            isp=gspecies(i,istr)
            if (isp.eq.1) then

                meam_f(i)=meame0(isp)*(rho_i(i)/meamrho0(isp))* &
                    log(rho_i(i)/meamrho0(isp))

            endif

        enddo

    endif

end subroutine embeddingfunction

    ! Code fragment to use if I code in look-up tables:
    !
    !      if (lookuptables.eqv..true.) then
    !
    !         print *,'using lookuptables for emb en calc'
    !
    !         do i=1,gn_inequivalentsites(istr)
    !            isp=gspecies(i,istr)
    !
    !            if (isp.eq.1) then
    !
    !              call splint(rho_embFe,embFe,secder_embFe,
    !     +              emb_Fe_arrsz,emb_Fe_arrsz,
    !     +              rho_i(i),meam_f(i)) !***need to put index in before
    !                                        !rij
    !
    !            elseif (isp.eq.2) then
    !
    !              call splint(rho_embC,embC,secder_embC,
    !     +              emb_C_arrsz,emb_C_arrsz,
    !     +              rho_i(i),meam_f(i)) !***need to put index in before
    !                                        !rij
    !
    !            else
    !               print *,'stop in emb en routine'
    !               stop
    !            endif
    !         enddo
    !
    !      else
subroutine extendlattice

    !----------------------------------------------------------------------c
    !
    !     Extend the lattice (I.e. repeat the supercell using the
    !     lattice vectors) so that when an atom is displaced to calculate
    !     an atomic force, the same atoms in neighboring supercells are not
    !     also displaced.
    !
    !     Called by: initializestruc
    !     Calls: calculateshells
    !     Returns: -
    !     Files read: -
    !     Files written: -
    !
    !     Andy Duff, Jan 2008
    !
    !----------------------------------------------------------------------c


    use m_generalinfo
    use m_datapoints
    use m_poscar
    use m_meamparameters
    use m_atomproperties

    implicit none
    logical inlist
    integer nshell,nnatoms_tot,i,it1,it2,it3,j,nnatoms
    integer, allocatable:: zz_tmp(:)
    real(8) rmax2,dd, &
        dxyz(3),xyztry1(3),xyztry2(3),xyztry(3)
    real(8), allocatable:: coords_tmp(:,:)

    !Store old coordinates and z values
    norigatoms=natoms
    allocate(coords_tmp(1:3,natoms),zz_tmp(natoms))
    do i=1,natoms
        coords_tmp(1:3,i)=coordinates(1:3,i)
        zz_tmp(i)=zz(i)
    enddo

    !Determine how many shells of unit cells need to be included
    !around the central unit cell in order that all atoms within
    !p_rmax of any atom in the central unit cell are included.
    call calculateshells(nshell,nnatoms_tot,p_rmax)

    !Set-up new coordinates
    deallocate(coordinates,zz)
    allocate(coordinates(3,nnatoms_tot),zz(nnatoms_tot))
    do i=1,natoms
        coordinates(1:3,i)=coords_tmp(1:3,i)
        zz(i)=zz_tmp(i)
    enddo
    deallocate(coords_tmp,zz_tmp)

    !For each atom in the additional shells of unit cells determined by
    !'calculateshells', see which atoms are within p_max of one of the
    !inequivalent atoms (I.e. those in the 'central' cell) - and only retain
    !those atoms.
    rmax2=p_rmax**2
    nnatoms=natoms
    do i=1,natoms
        do it1=-nshell,nshell
            xyztry1=coordinates(1:3,i)+tr(1:3,1)*real(it1)
            do it2=-nshell,nshell
                xyztry2=xyztry1+tr(1:3,2)*real(it2)
                do it3=-nshell,nshell
                    xyztry=xyztry2+tr(1:3,3)*real(it3)
                    inlist=((it1.eq.0).and.(it2.eq.0).and.(it3.eq.0))
                    do j=1,natoms
                        dxyz=xyztry-coordinates(1:3,j)
                        dd=dxyz(1)**2+dxyz(2)**2+dxyz(3)**2
                        if ((dd.le.rmax2).and.(.not.inlist)) then
                            nnatoms=nnatoms+1 !add to the list
                            coordinates(1:3,nnatoms)=xyztry(1:3)
                            zz(nnatoms)=zz(i)
                            nat(z2species(zz(i)))=nat(z2species(zz(i)))+1
                            inlist=.true.
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    natoms=nnatoms

end subroutine
subroutine extractlmax

    !--------------------------------------------------------------c
    !
    !     Read in lmax from a user supplied potential file.
    !
    !     Called by:     program MEAMfit
    !     Returns:       lmax
    !     Files read:    startparameterfile
    !
    !     Andrew Duff 2014
    !
    !--------------------------------------------------------------c

    use m_filenames
    use m_meamparameters
    use m_atomproperties

    implicit none

    open(unit=1,file=trim(startparameterfile),status='old')
    read(1,*) lmax
    close(1)

end subroutine extractlmax
subroutine findsmallestsepn

    !     Calculate the smallest and largest seperations (less than p_rmax)
    !     as a function of pairs of species. Write a file containing the 
    !     smallest interatomic separation as a function of file number and 
    !     species, and also store for later output in the 'fitted_quantities'
    !     file. This data is useful when choosing fitting database and/or 
    !     understanding why some energies/forces are not well reproduced by 
    !     potential.
    !
    !     Andrew Duff, Imperial-College, 2014

    use m_geometry
    use m_atomproperties
    use m_generalinfo

    implicit none

    integer i,j,jj,nni,isp,jsp,cnt, &
        isp_forwrite,jsp_forwrite,ibin,isp_test,jsp_test
    integer, parameter :: nbin = 100 ! Num bins used to store atom sepns
    integer storeids(maxspecies,maxspecies,nstruct,2), &
        numsepns(maxspecies,maxspecies),numsepns_linarray(10)
    real(8) dsepn,smlsepn_overall,lrgsepn_overall,sepn
    real(8) smallestsepn,largestsepn,smlsepnperspc(maxspecies,maxspecies), &
        lrgsepnperspc(maxspecies,maxspecies),tofileSml(10),tofileLrg(10), &
        tofile(10)

    if (allocated(lrgsepnperspc_overall)) deallocate(lrgsepnperspc_overall)
    if (allocated(smlsepnperspc_overall)) deallocate(smlsepnperspc_overall)
    if (allocated(sdSepn)) deallocate(sdSepn,avgSepn,nSepn)

    allocate(smlsepnperspc_overall(maxspecies,maxspecies), &
             lrgsepnperspc_overall(maxspecies,maxspecies))
    allocate(sdSepn(maxspecies,maxspecies), &
             avgSepn(maxspecies,maxspecies), &
             nSepn(maxspecies,maxspecies))

    !Record largest and smallest seperations for each structure (and also the
    !same, for each species combinations). Also record the largest and smallest
    !seperations over all structures.
    open(60,file='smallestsepnvsstruc.dat')
    write(60,*) '# structure-id | smallest sepn'

    !Set-up initial values of smallest and largest overall (across all
    !structures) seperations to ensure, e.g., zero is not wrongly returned as
    !the smallest seperation.
    smlsepn_overall=1000d0
    lrgsepn_overall=0d0
    do isp=1,maxspecies
        do jsp=1,maxspecies
            smlsepnperspc_overall(isp,jsp)=1000d0
            lrgsepnperspc_overall(isp,jsp)=0d0
        enddo
    enddo
    numsepns=0d0

    do istr=1,nstruct
        !Set up initial values as above but for use 'per structure'
        smallestsepn=1000d0
        largestsepn=0d0
        do isp=1,maxspecies
            do jsp=1,maxspecies
                smlsepnperspc(isp,jsp)=1000d0
                lrgsepnperspc(isp,jsp)=0d0
            enddo
        enddo
        do i=1,gn_inequivalentsites(istr)
            isp=gspecies(i,istr)
            nni=gn_neighbors(i,istr)
            do jj=1,nni
                j=gneighborlist(jj,i,istr)
                jsp=gspecies(j,istr)
                smallestsepn=min(smallestsepn,diststr(jj,i,0,0,istr))
                largestsepn=max(largestsepn,diststr(jj,i,0,0,istr))
                smlsepn_overall=min(smlsepn_overall, &
                    diststr(jj,i,0,0,istr))
                lrgsepn_overall=max(lrgsepn_overall, &
                    diststr(jj,i,0,0,istr))
                !Ensure we only store separations between different species
                !in the array smlsepnperspc(:,:) for the case where the
                !first index is smaller than the second.
                if (isp.gt.jsp) then
                    isp_forwrite=jsp
                    jsp_forwrite=isp
                else
                    isp_forwrite=isp
                    jsp_forwrite=jsp
                endif
                numsepns(isp_forwrite,jsp_forwrite)= &
                     numsepns(isp_forwrite,jsp_forwrite)+1
                if (diststr(jj,i,0,0,istr).lt. &
                    smlsepnperspc(isp_forwrite,jsp_forwrite)) then
                    smlsepnperspc(isp_forwrite,jsp_forwrite)= &
                        diststr(jj,i,0,0,istr)
                    storeids(isp_forwrite,jsp_forwrite,istr,1)=i
                    storeids(isp_forwrite,jsp_forwrite,istr,2)=j
                endif
                if ((diststr(jj,i,0,0,istr).gt. &
                    lrgsepnperspc(isp_forwrite,jsp_forwrite)) .and. &
                    (diststr(jj,i,0,0,istr).lt. &
                     p_rmax)) then
                    lrgsepnperspc(isp_forwrite,jsp_forwrite)= &
                        diststr(jj,i,0,0,istr)
                endif
                if (smallestsepn.eq.0d0) then
                    print *,'error, smallestsepn=0'
                    stop
                endif
            enddo
        enddo
        !Check if the smallest or largest seperations found for this structure
        !are the smallest or largest overall across all structures, and also
        !setup 'tofile' variables to help writing data to files
        cnt=0
        do isp=1,maxspecies
            do jsp=1,maxspecies
                if (isp.le.jsp) then
                    cnt=cnt+1
                    if (cnt.gt.10) then
                        print *,'increase array size in findsmallestsepn'
                        stop
                    endif
                    tofileSml(cnt)=smlsepnperspc(isp,jsp)
                    tofileLrg(cnt)=lrgsepnperspc(isp,jsp)
                    numsepns_linarray(cnt)=numsepns(isp_forwrite,jsp_forwrite)
                    smlsepnperspc_overall(isp,jsp)=min( &
                           smlsepnperspc_overall(isp,jsp), &
                           smlsepnperspc(isp,jsp) )
                    lrgsepnperspc_overall(isp,jsp)=max( &
                           lrgsepnperspc_overall(isp,jsp), &
                           lrgsepnperspc(isp,jsp) )
                endif
            enddo
        enddo
        smallestsepnStr(istr)=smallestsepn
        write(60,'(I4,A1,F10.5)',advance="no") istr,' ', &
            smallestsepn
        do i=1,cnt
            smlsepnperspcStr(i,istr)=tofileSml(i)
            lrgsepnperspcStr(i,istr)=tofileLrg(i)
            write(60,'(A1,F10.5)',advance="no") ' ',tofileSml(i)
        enddo
        smlsepnperspcStr_numEntries=cnt
        write(60,*)
    enddo

    write(60,*) 'smallest sepns across all files:'
    do isp=1,maxspecies
        do jsp=1,maxspecies
            if (isp.le.jsp) then
               write(60,*) isp,jsp,':',smlsepnperspc_overall(isp,jsp)
            endif
        enddo
    enddo
    write(60,*) 'overall :',smlsepn_overall
    write(60,*) 
    write(60,*) 'largest sepns across all files:'
    do isp=1,maxspecies
        do jsp=1,maxspecies
            if (isp.le.jsp) then
                write(60,*) isp,jsp,':',lrgsepnperspc_overall(isp,jsp)
            endif
        enddo
    enddo
    write(60,*) 'overall :',lrgsepn_overall
    close(60)

    !For each of the seperations provided in the previous file, record the id
    !numbers (in the order read in from the vasprun.xml/POSCAR files) for the
    !atoms corresponding to these seperations.
    open(60,file='idsforsmallestsepns')
    write(60,*)
    write(60,*) '# ids of atoms (i,j) for above smallest sepns'
    write(60,*) '(for species, e.g. (1,1); (1,2); (2,2))'
    do istr=1,nstruct
        write(60,'(I4,A2)',advance="no") istr,': '
        do isp=1,maxspecies
            do jsp=1,maxspecies
                if (isp.le.jsp) then
                    write(60,'(A1,I4,A1,I4)',advance="no") ' ', &
                        storeids(isp,jsp,istr,1),',',storeids(isp,jsp,istr,2)
                endif
            enddo
        enddo
        write(60,*)
    enddo
    close(60)

    !Also produce histogram showing distribution of separations for different
    !pairs of species
    open(60,file='sepnHistogram.out')
    write(60,*) '# sepn   no. sepns within range, per species (e.g. (1,1); (1,2); (2,2))'
    dsepn=(lrgsepn_overall-smlsepn_overall)/dble(nbin)
    do ibin=1,nbin
       sepn=smlsepn_overall+dble(ibin-1)*dsepn
       cnt=0
       do isp=1,maxspecies
          do jsp=1,maxspecies
             if (isp.le.jsp) then
                cnt=cnt+1
                tofile(cnt)=0
                do istr=1,nstruct
                   do i=1,gn_inequivalentsites(istr)
                       isp_test=gspecies(i,istr)
                       nni=gn_neighbors(i,istr)
                       do jj=1,nni
                           j=gneighborlist(jj,i,istr)
                           jsp_test=gspecies(j,istr)
                           if ((isp_test.eq.isp).and.(jsp_test.eq.jsp)) then
                              if ((diststr(jj,i,0,0,istr).gt.sepn).and. &
                                  (diststr(jj,i,0,0,istr).lt.sepn+dsepn)) then
                                 tofile(cnt)=tofile(cnt)+1
                              endif
                           endif
                       enddo
                   enddo
                enddo
             endif
          enddo
       enddo
       write(60,'(F10.5,A1)',advance="no") sepn,' '
       do i=1,cnt
           write(60,'(A1,F10.5)',advance="no") ' ',dble(tofile(i))/dble(numsepns_linarray(i))
       enddo
       write(60,*)
    enddo
    close(60)

    !Determine average seperations between species
    avgSepn=0
    nSepn=0
    do isp=1,maxspecies
       do jsp=1,maxspecies
          if (isp.le.jsp) then
             do istr=1,nstruct
                do i=1,gn_inequivalentsites(istr)
                   isp_test=gspecies(i,istr)
                   nni=gn_neighbors(i,istr)
                   do jj=1,nni
                      j=gneighborlist(jj,i,istr)
                      jsp_test=gspecies(j,istr)
                      if ((isp_test.eq.isp).and.(jsp_test.eq.jsp)) then
                         avgSepn(isp,jsp)=avgSepn(isp,jsp)+diststr(jj,i,0,0,istr)
                         nSepn(isp,jsp)=nSepn(isp,jsp)+1
                      endif
                   enddo
                enddo
             enddo
             avgSepn(isp,jsp)=avgSepn(isp,jsp)/nSepn(isp,jsp)
          endif
       enddo
    enddo

    sdSepn=0
    do isp=1,maxspecies
       do jsp=1,maxspecies
          if (isp.le.jsp) then
             do istr=1,nstruct
                do i=1,gn_inequivalentsites(istr)
                   isp_test=gspecies(i,istr)
                   nni=gn_neighbors(i,istr)
                   do jj=1,nni
                      j=gneighborlist(jj,i,istr)
                      jsp_test=gspecies(j,istr)
                      if ((isp_test.eq.isp).and.(jsp_test.eq.jsp)) then
                         sdSepn(isp,jsp)=sdSepn(isp,jsp)+ &
                   (diststr(jj,i,0,0,istr)-avgSepn(isp,jsp))**2
                      endif
                   enddo
                enddo
             enddo
             sdSepn(isp,jsp)=sqrt(sdSepn(isp,jsp)/nSepn(isp,jsp))
             !print *,'For isp=',isp,', jsp=',jsp,', average sepn=',avgSepn(isp,jsp), &
             !    ', s.d.=',sdSepn(isp,jsp)
          endif
       enddo
    enddo


end subroutine findsmallestsepn

subroutine Fnew(nfreep, popt, nf, F, uip, urp, ufp )

    !-------------------------------------------------------------------c
    !
    !     Calculate quality of the fit. The data-points are either
    !     meam energies or meam forces.
    !
    !     Called by:     optimizeparameters
    !     Calls:         parachange1, setupsplines, p_to_variables,
    !                 meamforce, meamenergy, getlargestbkgrnddens,
    !                 writemeamdatapoints, writederivedquantities,
    !                 dislaytime, plotfunctions
    !     Returns:       -
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Dec 2007-2015
    !
    !-------------------------------------------------------------------c

    use m_datapoints
    use m_geometry
    use m_optimization
    use m_generalinfo
    use m_meamparameters
    use m_filenames
    use m_poscar
    use m_electrondensity
    use m_atomproperties
    use m_plotfiles

    implicit none

    logical signflipped(np),minimumfailed
    character*80 filename,string1,string2
    integer i,j,ip,idatapoint,nfreep,foundM,ndatapointstofit,hrs
    integer nf,uip(*)
    real(8) funcforcdenom,funcendenom,F,Fen,Ffrc, &
        enbefore,enmiddle,enafter, &
        sumpp,summeamf
    !integer isp,jsp ! temp testing
    !real(8) sepn,pairpot ! temp variables for testing
    real(8) force(3),cutoff(np)
    real(8) urp(nfreep,3)
    real(8) popt(nfreep)

    external ufp

    !Convert parameters back to their correct values
    !(the popt array is used for the conjugate gradient optimizer, which prefers
    !numbers of a similar order)
    ip=0
    do i=1,np
        if ( (freep(i).eq.1).or.(freep(i).eq.2) ) then
            ip=ip+1
            if (p_orig(i).ne.0d0) then
                p(i)=popt(ip)*(p_orig(i)/poptStrt)
            else
                p(i)=popt(ip)
            endif
        endif
    enddo

    call parachange1(signflipped,cutoff,nf) !Ensure meam parameters have
                                            !the correct signs and are in the correct range
   !  call setupsplines !Prepare splines for pair-potentials and
   !                    !radial densities
    call p_to_variables !Convert variables in p() to 'sensibly-named' variables

   ! Temporary setup of array for testing purposes ----
   ! Turned out not to be faster than sums over cubic terms!!!
   ! do isp=1,maxspecies
   !   do jsp=1,maxspecies
   !      !print *,'isp=',isp,', jsp=',jsp
   !      if (jsp.ge.isp) then
   !         !print *,'i going form 1 to ',splnNvals
   !         do i=1,splnNvals
   !            sepn=cutoffMin+dble(i-1)*(cutoffMax-cutoffMin)/dble(splnNvals-1)
   !            pairpot=0d0
   !            do j=1,31,2
   !                if (sepn.le.pairpotparameter(j+1,isp,jsp)) then
   !                    pairpot=pairpot + pairpotparameter(j,isp,jsp) * &
   !                        ((pairpotparameter(j+1,isp,jsp)-sepn)**3)
   !                endif
   !            enddo
   !            pairpotStr(i,isp,jsp)=pairpot
   !            !if ((i.le.10).or.(i.ge.(splnNvals-10)).and.(isp.eq.1).and.(jsp.eq.1)) then
   !            !   write(81,*) sepn,pairpotStr(i,isp,jsp)
   !            !endif
   !         enddo
   !      endif
   !   enddo
   ! enddo
   ! stop
   ! print *,'finished spline init'
   ! --------------------------------------------------


    F=0d0
    Fen=0d0
    Ffrc=0d0
    funcforc=0d0
    funcen=0d0
    funcforcdenom=0d0
    funcendenom=0d0
    idatapoint=1
    foundM=0
    minimumfailed=.false.
    ndatapointstofit=0

    do i=1,nstruct

        if (optimizeforce(i).gt.0) then
            do j=1,gn_forces(i)

                !Compute 'funcforc', the s.d. of the force components - not
                !connected to the optimization function, computed further down
                if ((optforce(j,i).eqv..true.).and.(weights(i).ne.0d0)) &
                    then
                    !print *,'p(107)=',p(107),' p(108)=',p(108)
                    call meamforce(i,j,force) !Calculate the force acting
                    !on atom
                    fitdata(idatapoint:idatapoint+2)=force(1:3)
                    if (optimizeforce(i).eq.1) then
                        funcforc=funcforc+ &
                            ((fitdata(idatapoint)-truedata(idatapoint))**2 &
                            +(fitdata(idatapoint+1)-truedata(idatapoint+1))**2 &
                            +(fitdata(idatapoint+2)-truedata(idatapoint+2))**2 &
                            )
                    elseif (optimizeforce(i).eq.2) then
                        !Following code kept for legacy purposes, but not
                        !recommended to use:
                        if (truedata(idatapoint).ne.0d0) then
                            funcforc=funcforc+ &
                                ((fitdata(idatapoint)-truedata(idatapoint))/ &
                                truedata(idatapoint))**2
                        endif
                        if (truedata(idatapoint+1).ne.0d0) then
                            funcforc=funcforc+ &
                                ((fitdata(idatapoint+1)-truedata(idatapoint+1))/ &
                                truedata(idatapoint+1))**2
                        endif
                        if (truedata(idatapoint+2).ne.0d0) then
                            funcforc=funcforc+ &
                                ((fitdata(idatapoint+2)-truedata(idatapoint+2))/ &
                                truedata(idatapoint+2))**2
                        endif
                    endif
                    funcforcdenom=funcforcdenom+3
                    ndatapointstofit=ndatapointstofit+3
                else
                    force=0d0
                    fitdata(idatapoint:idatapoint+2)=force(1:3)
                endif

                if (optforce(j,i).eqv..true.) then
                    !Add contributions to the optimization function due to forces. Here
                    !add to Ffrc, which will later be added to F.
                    if (optimizeforce(i).eq.1) then
                        F=F+weights(i)* &
                            ((fitdata(idatapoint)-truedata(idatapoint))**2 &
                            +(fitdata(idatapoint+1)-truedata(idatapoint+1))**2 &
                            +(fitdata(idatapoint+2)-truedata(idatapoint+2))**2 &
                            )
                        Ffrc=Ffrc+weights(i)* &
                            ((fitdata(idatapoint)-truedata(idatapoint))**2 &
                            +(fitdata(idatapoint+1)-truedata(idatapoint+1))**2 &
                            +(fitdata(idatapoint+2)-truedata(idatapoint+2))**2 &
                            )
                    elseif (optimizeforce(i).eq.2) then
                        !Following code kept for legacy purposes, but not
                        !recommended to use:
                        print *,'ERROR: fitting to relative errors of forces currently'
                        print *,'not supported, STOPPING.'
                        stop
                        if (truedata(idatapoint).ne.0d0) then
                            F=F+weights(i)* &
                                ((fitdata(idatapoint)-truedata(idatapoint))/ &
                                truedata(idatapoint))**2
                        endif
                        if (truedata(idatapoint+1).ne.0d0) then
                            F=F+weights(i)* &
                                ((fitdata(idatapoint+1)-truedata(idatapoint+1))/ &
                                truedata(idatapoint+1))**2
                        endif
                        if (truedata(idatapoint+2).ne.0d0) then
                            F=F+weights(i)* &
                                ((fitdata(idatapoint+2)-truedata(idatapoint+2))/ &
                                truedata(idatapoint+2))**2
                        endif
                    endif

                endif
                idatapoint=idatapoint+3

            enddo
            !            funcforcdenom=funcforcdenom+gn_forces(i)
        else

            if (ensureminimum(i).eqv..true.) then
                foundM=foundM+1
            endif

            if (ensureminimum(i).eqv..true.) then
                print *,'ERROR: ensureminimum mode currently ot supported,'
                print *,'STOPPING.'
                stop
                if (foundM.eq.1) then
                    call meamenergy(i,0,0,enbefore,sumpp,summeamf)
                elseif (foundM.eq.2) then
                    call meamenergy(i,0,0,enmiddle,sumpp,summeamf)
                elseif (foundM.eq.3) then
                    call meamenergy(i,0,0,enafter,sumpp,summeamf)
                endif
                fitdata(idatapoint)=0d0
            else
                if (weights(i).ne.0d0) then
                    call meamenergy(i,0,0,fitdata(idatapoint),sumpp,summeamf)
                    F=F+weights(i)* &
                        (fitdata(idatapoint)-truedata(idatapoint))**2
                    !Add energy contributions to optimization function (Fen,
                    !later to be added to F), as well as the s.d. of the
                    !energies (funcen)
                    if (useRef) then
                       Fen=Fen+weights(i)* &
                           ( (fitdata(idatapoint)-fitdata(1)) - (truedata(idatapoint)-truedata(1)) )**2
                       funcen=funcen+ &
                           ( (fitdata(idatapoint)-fitdata(1)) - (truedata(idatapoint)-truedata(1)) )**2
                    else
                       Fen=Fen+weights(i)* &
                           ( fitdata(idatapoint) - truedata(idatapoint) )**2
                       funcen=funcen+ &
                           ( fitdata(idatapoint) - truedata(idatapoint) )**2
                    endif
                    funcendenom=funcendenom+1
                    ndatapointstofit=ndatapointstofit+1
                else
                    fitdata(idatapoint)=0d0
                endif
            endif

            if ( (ensureminimum(i).eqv..true.).and. &
                (foundM.eq.3) ) then
                !Legacy code, not currently used:
                foundM=0
                if (enafter.lt.enmiddle) then
                    F=F+weights(i-2)*((enmiddle-enafter)**2)
                endif
                if (enbefore.lt.enmiddle) then
                    F=F+weights(i-2)*((enmiddle-enbefore)**2)
                endif
            endif
            idatapoint=idatapoint+1

        endif

    enddo

    i=i-1

   !!Old definition of optimization function
   !if ((F.lt.0d0).or.(F.ge.0d0)) then
   !    F=sqrt(F/real(ndatapointstofit))
   !else
   !    F=1000000000000d0
   !    if (minimumfailed.eqv..true.) then
   !        print *,'minimum condition failed'
   !    endif
   !endif

    !New definition of optimization function
    if ( ((Fen.lt.0d0).or.(Fen.ge.0d0)).and. &
         ((Ffrc.lt.0d0).or.(Ffrc.ge.0d0)) ) then
        F=0d0
        if (Fen.gt.0d0) then
           F=F+sqrt(Fen/(real(nEn)*varEn))
        endif
        if (Ffrc.gt.0d0) then
           F=F+sqrt(Ffrc/(real(nFrcComp)*varFrc))
           !print *,'Ffrc=',Ffrc,', nFrcComp=',nFrcComp,', varFrc=',varfrc,', Ffrc/(nFrcComp*varFrc)=',Ffrc/(real(nFrcComp)*varFrc),', F=',F
           !print *,'F=',F
        endif
     !   !New penalty for cut-offs exceeding bounds:
     !   call CalcCutoffPenalty
     !  if (CutoffPenalty.ne.0d0) then
     !   print *,'CutoffPenalty=',CutoffPenalty,' (compared to F W/O = ',F,')'
     !  stop
     !  endif
     !   F=F+CutoffPenalty
    else
        F=1000000000000d0
        if (minimumfailed.eqv..true.) then
            print *,'minimum condition failed'
        endif
    endif

    !Divide through funcforc and funcen to get s.d.s
    if (funcforcdenom.ne.0d0) then
        funcforc=sqrt(funcforc/funcforcdenom)
    endif
    if (funcendenom.ne.0d0) then
        funcen=sqrt(funcen/funcendenom)
    endif

    if (nsteps.eq.1) then

        !Since code is stopping now, write out all energies and forces
        !(including those which were given zero weight)
        idatapoint=1
        !Record the largest background density encountered (this
        !will inform the maxiumum density used when outputting the
        !lammps input file).
        open(60,file='largestrho.dat')
        write(60,*) '# structure-id | largest rho'
        maxlargestrho=0d0
        do i=1,nstruct
            if (optimizeforce(i).gt.0) then
                do j=1,gn_forces(i)
                    if (optforce(j,i).eqv..true.) &
                        then
                        call meamforce(i,j,force)
                        fitdata(idatapoint:idatapoint+2)=force(1:3)
                        !                  print *,idatapoint,force(1:3)
                    endif
                    idatapoint=idatapoint+3
                enddo
            else
                if (ensureminimum(i).eqv..false.) then
                    call meamenergy(i,0,0,fitdata(idatapoint), &
                                    sumppStr(idatapoint),summeamfStr(idatapoint))
                endif
                idatapoint=idatapoint+1
            endif
            call getlargestbkgrnddens
            maxlargestrho=MAX(maxlargestrho,largestrho)
            write(60,*) i,largestrho
        enddo
        print *
        print *,'-----------------------------------------'
        print *,'Optimization function=',F
        print *,'-----------------------------------------'
        print *
        print *,'rms error on energies=',funcen
        print *,'rms error on forces=',funcforc
        call writemeamdatapoints(6,.false.)
        !call writederivedquantities(6)
        open(61,file='fitted_quantities.out')
        call writemeamdatapoints(61,.true.)
        close(61)

        110     format (F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2, &
            A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2, &
            A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2, &
            A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2, &
            A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2, &
            A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2)
        close(51)

        print *
        print *,'----------------------'
        print *,'Optimization completed'
        call displaytime(6,hrs)
        print *,'----------------------'

        if (lmax.eq.0) then
           !Open file for output of potential in LAMMPS format
           filename=""
           do i=1,maxspecies
               write(string1,'(A2)') element(speciestoZ(i))
               filename=trim(filename)//trim(string1)
           enddo
           string2='.eam.alloy'
           filename=trim(filename)//trim(string2)
           open(58,file=trim(adjustl(filename)))
        endif

        !Output files (inc. LAMMPS, Camelon and potential plots
        call plotfunctions(.true.)

        if (lmax.eq.0) then
           close(58)
        endif

        if (verbose) then
           print *,'maximum background density encountered:',maxlargestrho
           print *,'(recorded per structure in largestrho.dat)'
        endif
        print *

        stop

    endif

end subroutine Fnew
subroutine generaterandomstartparas

    !--------------------------------------------------------------c
    !
    !     Randomly initializes the MEAM parameters according to the
    !     limits previously determined.
    !
    !     Called by: MEAMfit
    !     Calls: RANDOM_NUMBER
    !     Returns: cmin, cmax, meamtau, meamrhodecay, meame0,
    !         meamrho0, meamemb3, meamemb4, pairpotparameter,
    !         enconst
    !     Files read: -
    !     Files written: -
    !
    !     Andy Duff, Feb 2013.
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_filenames
    use m_atomproperties
    use m_generalinfo
    use m_optimization

    implicit none

    !logical sorted
    !integer iseed
    integer i,j,ii,ll,ispc,jspc,spc_cnt
    real(8) tmp,maxradius_current

    !Variables just used for reading in start_meam_parameters file:
    integer nfreep

    !---- cmin and cmax ----
    j=1
    ispc=1
    jspc=1
    do i=1,m3
        if (freep(i).ge.2) then
            if (freep(i+m3).lt.2) then
                do
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    cmin(j,ispc,jspc)=5d0*tmp
                    if (cmin(j,ispc,jspc).lt.cmax(j,ispc,jspc)) then
                        exit
                    endif
                enddo
            else
                !tmp=ran(iseed)
                call RANDOM_NUMBER(tmp)
                cmin(j,ispc,jspc)=5d0*tmp
            endif
        endif
        if (j.lt.m1) then
            j=j+1
        else
            if (ispc.lt.m1) then
                j=1
                ispc=ispc+1
            else
                j=1
                ispc=1
                jspc=jspc+1
            endif
        endif
    enddo
    j=1
    ispc=1
    jspc=1
    do i=m3+1,2*m3
        if (freep(i).ge.2) then
            do
                !tmp=ran(iseed)
                call RANDOM_NUMBER(tmp)
                cmax(j,ispc,jspc)=5d0*tmp
                if (cmin(j,ispc,jspc).lt.cmax(j,ispc,jspc)) then
                    exit
                endif
            enddo
        endif
        if (j.lt.m1) then
            j=j+1
        else
            if (ispc.lt.m1) then
                j=1
                ispc=ispc+1
            else
                j=1
                ispc=1
                jspc=jspc+1
            endif
        endif
    enddo
    !----------------------

    !---- meamtau ----
    if (lmax.gt.0) then
       do ll=1,lmax
           do ispc=1,maxspecies
               !Check if we need to generate random paramaters for
               !meamtau(ll,isp)
               if (freep(2*m3+ispc+ll*m1).ge.2) then
                   !tmp=ran(iseed)
                   call RANDOM_NUMBER(tmp)
                   !print *,'rand nmr=',tmp
                   meamtau(ll,ispc)=10d0**(meamtau_minorder(ll,ispc)+tmp* &
                       (meamtau_maxorder(ll,ispc)-meamtau_minorder(ll,ispc)))
                   !tmp=ran(iseed)
                   call RANDOM_NUMBER(tmp)
                   !print *,'rand nmr=',tmp
                   if (tmp.gt.0.5d0) then
                       meamtau(ll,ispc)=-meamtau(ll,ispc)
                   endif
                   !print *,'meamtau(',ispc,',',ll,')=',meamtau(ll,ispc)
               endif
           enddo
       enddo
    endif
    !-----------------

    !---- electron density functions ----
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    !Check if we need to generate random paramaters for
                    !thi(ispc,jspc,ll)
                    nfreep=0
                    do i=1,12
                        if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i).ge.2) then
                            nfreep=nfreep+1
                        endif
                    enddo
                    if (nfreep.gt.0) then
                        if (typethi(ispc,jspc).eq.1) then !Cubic form

                            !Generate randomly generated parameters.
                            !Define an overall maximum, maxradius_current, for which
                            !all cut-off radii must be less than (avoids
                            !thi(l,ispc,jspc) preferentially initializing with a large
                            !radial extent.
                            !tmp=ran(iseed)       !meamrhodecay_minradius(ll,ispc,jspc)
                            call RANDOM_NUMBER(tmp)
                            !print *,'rand nmr=',tmp
                            maxradius_current=meamrhodecay_minradius(ll,ispc,jspc)+ &
                                0.1d0*(meamrhodecay_maxradius(ll,ispc,jspc)- &
                                       meamrhodecay_minradius(ll,ispc,jspc)) + &
                                      (meamrhodecay_maxradius(ll,ispc,jspc)- &
                                       (meamrhodecay_minradius(ll,ispc,jspc)+ &
                                        0.1d0*(meamrhodecay_maxradius(ll,ispc,jspc)- &
                                        meamrhodecay_minradius(ll,ispc,jspc))))*tmp
                            do i=1,6
                                if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+2*i).ge.2) &
                                    then
                                    !tmp=ran(iseed)
                                    call RANDOM_NUMBER(tmp)
                                    !print *,'rand nmr=',tmp
                                    meamrhodecay(2*i,ll,ispc,jspc)= &
                                       meamrhodecay_minradius(ll,ispc,jspc)+ &
                                  !     (maxradius_current- &
                                        (meamrhodecay_maxradius(ll,ispc,jspc)- &
                                        meamrhodecay_minradius(ll,ispc,jspc))*tmp
                                    !print *,'meamrhodecay(',2*i,',',ll,',',ispc,',',jspc,')=',meamrhodecay(2*i,ll,ispc,jspc)
                                endif
                                !Put -1d0 placeholders into the coefficient slots that
                                !need to be randomly generated, then after re-ordering
                                !we know which coefficients to generate.
                                if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+2*i-1).ge.2) &
                                    then
                                    meamrhodecay(2*i-1,ll,ispc,jspc)=-1d0
                                endif
                            enddo
                          ! print *,'random cutoffs for elec_Zr(',ispc,',',jspc,'):',meamrhodecay(2,0,ispc,jspc),meamrhodecay(4,0,ispc,jspc)
                          ! stop

                            !Generate random coefficients
                            ii=1
                            do i=1,6
                                if (meamrhodecay(2*i-1,ll,ispc,jspc).eq.-1d0) then
                                    !tmp=ran(iseed)
                                    call RANDOM_NUMBER(tmp)
                                    !print *,'rand nmr=',tmp
                                    meamrhodecay(2*i-1,ll,ispc,jspc)= &
                                        10d0**(meamrhodecay_minorder(ii,ll,ispc,jspc)+ &
                                               tmp*(meamrhodecay_maxorder(ii,ll,ispc,jspc)- &
                                                    (meamrhodecay_minorder(ii,ll,ispc,jspc))))
                                    if (meamrhodecay_negvals(ll,ispc,jspc).eq.2) then
                                        meamrhodecay(2*i-1,ll,ispc,jspc)= &
                                            -meamrhodecay(2*i-1,ll,ispc,jspc)
                                    elseif (meamrhodecay_negvals(ll,ispc,jspc).eq.3) then
                                        !tmp=ran(iseed)
                                        call RANDOM_NUMBER(tmp)
                                        if (tmp.gt.0.5d0) then
                                            meamrhodecay(2*i-1,ll,ispc,jspc)= &
                                                -meamrhodecay(2*i-1,ll,ispc,jspc)
                                        endif
                                    endif
                                    !print *,'meamrhodecay(',2*i-1,',',ll,',',ispc,',',jspc,')=',meamrhodecay(2*i-1,ll,ispc,jspc)
                                    ii=ii+1
                                endif
                            enddo
                        else
                            print *,'     typethi(',ispc,',',jspc,')=', &
                                typethi(ispc,jspc),'not supported, stopping'
                            stop
                        endif
                    endif
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo
    !------------------------------------


    !---- embedding functions ----
    do ispc=1,maxspecies
        do i=1,4
            !print *,'For the embedding functions, freep(',2*m3+lm1*m1+12*lm1*m2+ispc+(i-1)*m1,', which corresponds to ispc=',ispc,' and i=',i,' has the value:',freep(2*m3+lm1*m1+12*lm1*m2+ispc+(i-1)*m1),' (with embfuncRandType=',embfuncRandType,')'
            if (freep(2*m3+lm1*m1+12*lm1*m2+ispc+(i-1)*m1).ge.2) then
               if (embfuncRandType.eq.1) then
                if (i.eq.1) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    !print *,'rand nmr=',tmp
                    meame0(ispc)=tmp*(meame0_maxorder(ispc)- &
                        meame0_minorder(ispc))+meame0_minorder(ispc)
                    !print *,'meame0(',ispc,')=',meame0(ispc)
                elseif (i.eq.2) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    !print *,'rand nmr=',tmp
                    meamrho0(ispc)=tmp*(meamrho0_maxorder(ispc)- &
                        meamrho0_minorder(ispc))+meamrho0_minorder(ispc)
                    !print *,'meamrho0(',ispc,')=',meamrho0(ispc)
                elseif (i.eq.3) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    !print *,'rand nmr=',tmp
                    meamemb3(ispc)=tmp*(meamemb3_maxorder(ispc)- &
                        meamemb3_minorder(ispc))+meamemb3_minorder(ispc)
                    !print *,'meamemb3(',ispc,')=',meamemb3(ispc)
                elseif (i.eq.4) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    meamemb4(ispc)=tmp*(meamemb4_maxorder(ispc)- &
                        meamemb4_minorder(ispc))+meamemb4_minorder(ispc)
                endif
               elseif (embfuncRandType.eq.2) then
                if (i.eq.1) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    !print *,'rand nmr=',tmp
                    meame0(ispc)=10d0**( tmp*(meame0_maxorder(ispc)- &
                        meame0_minorder(ispc))+meame0_minorder(ispc) )
                    !print *,'meame0(',ispc,')=',meame0(ispc)
                elseif (i.eq.2) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    !print *,'rand nmr=',tmp
                    meamrho0(ispc)=10d0**( tmp*(meamrho0_maxorder(ispc)- &
                        meamrho0_minorder(ispc))+meamrho0_minorder(ispc) )
                    !print *,'meamerho0(',ispc,')=',meame0(ispc)
                elseif (i.eq.3) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    !print *,'rand nmr=',tmp
                    meamemb3(ispc)=10d0**( tmp*(meamemb3_maxorder(ispc)- &
                        meamemb3_minorder(ispc))+meamemb3_minorder(ispc) )
                    !print *,'meamemb3(',ispc,')=',meamemb3(ispc)
                elseif (i.eq.4) then
                    !tmp=ran(iseed)
                    call RANDOM_NUMBER(tmp)
                    meamemb4(ispc)=10d0**( tmp*(meamemb4_maxorder(ispc)- &
                        meamemb4_minorder(ispc))+meamemb4_minorder(ispc) )
                endif
               else
                print *,'ERROR: method for randomly generating embfuncs not supported (embfuncRandType=',embfuncRandType,')'
               endif
            endif
        enddo
    enddo
    !----------------------------

    !---- pair-potentials ----
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                if (nfreeppairpot(ispc,jspc).gt.0) then
                    if (typepairpot(ispc,jspc).eq.2) then
                        !Generate randomly generated parameters.
                        !Define an overall maximum, maxradius_current, for which
                        !all cut-off radii must be less than (avoids
                        !thi(l,ispc,jspc) preferentially initializing with a
                        !large radial extent.
                        !tmp=ran(iseed)
                        call RANDOM_NUMBER(tmp)
                        !print *,'rand nmr=',tmp
                        maxradius_current=pairpotparameter_minradius(ispc,jspc)+ &
                            0.1d0*(pairpotparameter_maxradius(ispc,jspc)- &
                                   pairpotparameter_minradius(ispc,jspc)) + &
                            (pairpotparameter_maxradius(ispc,jspc)- &
                             (pairpotparameter_minradius(ispc,jspc)+0.1d0* &
                              (pairpotparameter_maxradius(ispc,jspc)- &
                               pairpotparameter_minradius(ispc,jspc)) ) )*tmp
                        do i=1,16
                            if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+2*i+32*spc_cnt) &
                                .ge.2) then
                                !tmp=ran(iseed)
                                call RANDOM_NUMBER(tmp)
                                !print *,'rand nmr=',tmp
                                pairpotparameter(2*i,ispc,jspc)= &
                                    pairpotparameter_minradius(ispc,jspc)+ &
                            !        (maxradius_current-pairpotparameter_minradius(ispc,jspc))*tmp
                                    (pairpotparameter_maxradius(ispc,jspc)-pairpotparameter_minradius(ispc,jspc))*tmp
                            ! print *,'random cutoffs for pairpot cutoff(',2*i,',',ispc,',',jspc,'):',pairpotparameter(2*i,ispc,jspc)
                                !print *,'pairpotparameter(',2*i,',',ispc,',',jspc,')=',pairpotparameter(2*i,ispc,jspc)
                            endif

                            !Put -1d0 placeholders into the coefficient slots that
                            !need
                            !to be randomly generated. This is so that, after
                            !re-ordering, we know which coefficients to generate.
                            if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+2*i-1+32*spc_cnt) &
                                .ge.2) then
                                pairpotparameter(2*i-1,ispc,jspc)=-1d0
                            endif
                        enddo
                        !Sort the terms in order of increasing cut-off radii
                !       do
                !           sorted=.true.
                !           do i=2,30,2
                !               if ((pairpotparameter(i+2,ispc,jspc).lt. &
                !                   pairpotparameter(i,ispc,jspc)).and. &
                !                   (pairpotparameter(i+2,ispc,jspc).ne.0d0)) then
                !                   tmp=pairpotparameter(i-1,ispc,jspc)
                !                   tmp2=pairpotparameter(i,ispc,jspc)
                !                   pairpotparameter(i-1,ispc,jspc)= &
                !                       pairpotparameter(i+1,ispc,jspc)
                !                   pairpotparameter(i,ispc,jspc)= &
                !                       pairpotparameter(i+2,ispc,jspc)
                !                   pairpotparameter(i+1,ispc,jspc)=tmp
                !                   pairpotparameter(i+2,ispc,jspc)=tmp2
                !                   sorted=.false.
                !               endif
                !           enddo
                !           if (sorted.eqv..true.) then
                !               exit
                !           endif
                !       enddo
                        !Generate random coefficients
                        ii=1
                        do i=1,16
                            if (pairpotparameter(2*i-1,ispc,jspc).eq.-1d0) then
                                !tmp=ran(iseed)
                                call RANDOM_NUMBER(tmp)
                                !print *,'rand nmr=',tmp
                                pairpotparameter(2*i-1,ispc,jspc)= &
                                    10d0**(pairpotparameter_minorder(ii,ispc,jspc)+ &
                                    tmp*(pairpotparameter_maxorder(ii,ispc,jspc)- &
                                    (pairpotparameter_minorder(ii,ispc,jspc))))
                                if (pairpotparameter_negvals(ispc,jspc).eq.2) then
                                    pairpotparameter(2*i-1,ispc,jspc)= &
                                        -pairpotparameter(2*i-1,ispc,jspc)
                                elseif (pairpotparameter_negvals(ispc,jspc).eq.3) then
                                    !tmp=ran(iseed)
                                    call RANDOM_NUMBER(tmp)
                                    !print *,'rand nmr=',tmp
                                    if (tmp.gt.0.5d0) then
                                        pairpotparameter(2*i-1,ispc,jspc)= &
                                            -pairpotparameter(2*i-1,ispc,jspc)
                                    endif
                                endif
                                !print *,'pairpotparameter(',2*i-1,',',ispc,',',jspc,')=',pairpotparameter(2*i-1,ispc,jspc)
                                ii=ii+1
                            endif
                        enddo
                    else
                        print *,'typepairpot=',typepairpot(ispc,jspc), &
                            'not supported,stopping'
                        stop
                    endif
                endif
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo
    !--------------------------------------------------------
!    stop

    !---- enconst ----
    do ispc=1,maxspecies
        if (freep(m4+2*m2+ispc).ge.2) then
            !tmp=ran(iseed)
            call RANDOM_NUMBER(tmp)
            !print *,'rand nmr=',tmp
            enconst(ispc)=10d0**( tmp*(enconst_maxorder(ispc)-enconst_minorder(ispc)) + &
                               enconst_minorder(ispc) )
            !tmp=ran(iseed)
            call RANDOM_NUMBER(tmp)
            !print *,'rand nmr=',tmp
            if (tmp.gt.0.5d0) then
                enconst(ispc)=-enconst(ispc)
            endif
            !print *,'enconst(',ispc,')=',enconst(ispc)
        endif
    enddo
    !------------------

end subroutine generaterandomstartparas
subroutine getlargestbkgrnddens

    !----------------------------------------------------------------------c
    !
    !     Determines the largest value of background density (rho)
    !     encountered across all structures.
    !
    !     Called by:     optimizeparameters, fit_quality
    !     Returns:       largestrho
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Oct 2014
    !
    !----------------------------------------------------------------------c
 

    use m_electrondensity
    use m_geometry
    use m_plotfiles

    implicit none

    integer i

    largestrho=0d0
    do i=1,gn_inequivalentsites(istr)
        largestrho=MAX(largestrho,rho_i(i))
    enddo

end subroutine getlargestbkgrnddens
subroutine getnumconfigs(string,limits,nlimits,step)

    !-------------------------------------------------------------------c
    !
    !     Extract from 'string' the ranges of ionic configurations which 
    !     are to be used in the optimization, and then assign the limits
    !     of these ranges to 'limits', with the number of ranges
    !     assigned to 'nlimits'. Also check for 's', I.e., specification
    !     of a 'step' to be taken between adjacent configurations.
    !     Return these in 'step'.
    !
    !     Called by:     setupfilenames
    !     Calls:         checkifnumber
    !     Returns:       limits, nlimits, step
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff 2014
    !
    !-------------------------------------------------------------------c

    use m_optimization

    implicit none

    logical num,awaiting_upper,awaiting_step
    integer i,ilower,istruc,nlimits
    integer limits(10),step(10)
    character*1 char
    character*80 string,tmpstring

    i=1
    nlimits=0
    ilower=1
    awaiting_upper=.false.
    awaiting_step=.false.
    step=1 !Steps between adjancent configs set to 1 by default

    !Step through 'string' sequentially, checking if each character is a number
    !of a letter.
    do
       char=string(i:i)
       call checkifnumber(char,num)
       if (num.eqv..false.) then
          if (awaiting_step) then
             !If an 's' has previously been found, then the value of the step
             !has now been found
             awaiting_step=.false.
             tmpstring=trim(string(ilower:i-1))
             read(tmpstring,'(I10)') step(nlimits-1)
          else
             !...otherwise the presence of a character indicates that a limit
             !has been found (a lower or upper limit on configuration number)
             nlimits=nlimits+1
             if (nlimits.gt.10) then
                print *,'ERROR: more than 10 limits read in from fitdbse,'
                print *,'STOPPING'
                stop
             endif
             tmpstring=trim(string(ilower:i-1))
             read(tmpstring,'(I10)') limits(nlimits)
             if (char.eq.'-') then
                !We have found the lower limit...
                awaiting_upper=.true.
             elseif ((char.eq.';').or.(char.eq.' ').or.(char.eq.'s')) then
                if (awaiting_upper.eqv..false.) then
                   !We have found a single configuration
                   nlimits=nlimits+1
                   limits(nlimits)=limits(nlimits-1)
                   if (char.eq.'s') then
                      !Can't specify a step if we only have one configuration!
                      print *,'ERROR: step specified for single config in fitdbse,'
                      print *,'STOPPING.'
                      stop
                   endif
                else
                   !We have found the upper limit
                   awaiting_upper=.false.
                   if (char.eq.'s') then
                      awaiting_step=.true.
                   endif
                endif
             endif
          endif
          !Store the location for the start of the next string of numbers:
          ilower=i+1
       endif
       if (char.eq." ") then
          exit
       endif
       i=i+1
    enddo

end subroutine getnumconfigs
subroutine getRmax

    !----------------------------------------------------------------------c
    !
    !     Set-up the maximum radius to be considered subsequently. If
    !     cutoffs are to be optimized, this will be set either to 
    !     CUTOFFMAX, or if a potential file is provided by the user and
    !     some of the (fixed) cutoffs in this file are larger, then it will
    !     be set to this instead. If cutoffs are not to be optimized,
    !     it will either be set to the maximum cutoff in the potential
    !     input file, or if none is provided, it will be set to CUTOFFMAX
    !     (for the purposes of plotting the separation histogram).
    !
    !     Called by:     MEAMfit
    !     Calls:         -
    !     Returns:       p_rmax
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff, 2014
    !
    !----------------------------------------------------------------------c

    use m_meamparameters
    use m_generalinfo
    use m_atomproperties
    use m_optimization
    use m_filenames

    implicit none

    integer i,j,k,ispc,jspc,ll,spc_cnt
    real(8) rlargest

    if (readpotfile.eqv..false.) then
       !If no input potential file has been supplied the maximum radius to be
       !considered is just that given by CUTOFF_MAX
       if (verbose) then
          print *,'No potential parameters file supplied, therefore setting'
          print *,'p_rmax=CUTOFF_MAX (=',cutoffMax,')'
       endif
       p_rmax=cutoffMax
    endif

    !Check to see if any of the cut-off radii are being optimized, and keep
    !track of the largest of these.
    !First check the electron densities:
    cutoffopt=.false.
    spc_cnt=0
    rlargest=0d0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    do i=2,12,2
                        if ( (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i).eq.1).or. &
                             (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i).eq.2) ) then
                            cutoffopt=.true.
                        endif
                        if (meamrhodecay(i,ll,ispc,jspc).gt.rlargest) then
                           rlargest=meamrhodecay(i,ll,ispc,jspc)
                        endif
                    enddo
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo
    !Next the pair-potentials:
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                do i=2,32,2
                    if ( (freep(2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt).eq.1).or. &
                         (freep(2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt).eq.2) ) then
                        cutoffopt=.true.
                    endif
                    if (pairpotparameter(i,ispc,jspc).gt.rlargest) then
                       rlargest=pairpotparameter(i,ispc,jspc)
                    endif
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    if (verbose) then
       print *,'cutoffopt=',cutoffopt,', so...'
    endif
    if ((cutoffopt.eqv..false.).and.(readpotfile.eqv..true.)) then
       !In this case p_rmax can be set equal to the maximum radius used in the
       !start_meam_parameters file
       if (verbose) then
          print *,'Largest cutoff in start_meam_parameters=',rlargest
          print *,'Since these are not being optimized, setting p_rmax=',rlargest
       endif
       p_rmax=rlargest
       cutoffMax=rlargest !Even though we are not optimizing, this needs to be
                 !set to the largest cut-off to avoid the cut-offs being changed
    else
       if (cutoffMax.gt.rlargest) then
          if (verbose) then
             print *,'cutoffMax > largest cut-off in start_meam_parameters, therefore'
             print *,'setting p_rmax=',cutoffMax
          endif
          p_rmax=cutoffMax
       else
          if (verbose) then
             print *,'cutoffMax < largest cut-off in start_meam_parameters, therefore'
             print *,'setting p_rmax=',rlargest
          endif
          p_rmax=rlargest
       endif
    endif

end subroutine getRmax

subroutine getsettingspara(VarName,VarValue,foundvariable)

    !-------------------------------------------------------------------c
    !
    !     Read in the settings parameter specified by 'VarName'. If
    !     found, return foundvariable=.true.
    !
    !     Called by:     readsettings
    !     Calls:         optionalarg
    !     Returns:       VarValue, foundvariable
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff 2014
    !
    !-------------------------------------------------------------------c

    use m_filenames

    implicit none

    logical typelog,foundvariable,equals
    integer i,typeint,IOstatus,nwords,VarNamePlusEqLength,StringLength,wordLength
    real(8) typedble

    character*80 string
    character*80 VarName,word,VarValue,VarNamePlusEq
    character*80, allocatable:: words(:)
    character*1 tmpstring

    open(unit=1,file=trim(settingsfile))

    foundvariable=.false.
    VarNamePlusEq=trim(VarName)//"="
    do
       read(1,'(a80)',IOSTAT=IOstatus) string

       if (IOstatus.lt.0) then
          exit
       endif

       call optionalarg(string,nwords)
       if (nwords.gt.0) then

        allocate(words(1:nwords))
        read(string,*) words
        !print *,'string=',string

        !First check that line is not commented:
        if (words(1).ne.'#') then

          if (nwords.ge.1) then
             !Check if first word has an equals sign in it (eg,
             !VERBOSE=TRUE)
             word=words(1)
             wordLength=len_trim(word)
             equals=.false.
             do i=1,wordLength
               if (word(i:i).eq."=") then
                  equals=.true.
                  !print *,'equals in first word.'
                  exit
               endif
             enddo
             if (equals) then
               !print *,'checking if ',trim(word(1:i-1)), &
               !'equals',trim(VarName)
               if (trim(word(1:i-1)).eq.trim(VarName)) then
                 !print *,'i=',i,' wordLength=',wordLength
                 if (i.lt.wordLength) then
                   !Value of variable found in first word. Extract it:
                   !print *,'..found it within first word..'
                   VarValue=word(i+1:wordLength)
                   foundvariable=.true.
                   exit
                 else
                   if (nwords.ge.2) then
                     !Value of variable found in second word. Extract it:
                     !print *,'..found it in second word..'
                     VarValue=words(2)
                     foundvariable=.true.
                     exit
                   endif
                 endif
               endif
             else
               !print *,'checking if ',trim(words(1)), &
               !'equals',trim(VarName)
               if (trim(words(1)).eq.trim(VarName)) then
                 if (nwords.ge.2) then
                   if (trim(words(2)).eq."=") then
                     if (nwords.ge.3) then
                       !Value of variable found in third word. Extract it:
                       !print *,'..found it in third word..'
                       VarValue=words(3)
                       foundvariable=.true.
                       exit
                     endif
                   else
                     word=words(2)
                     wordLength=len_trim(word)
                     do i=1,wordLength
                       if (word(i:i).eq."=") then
                         equals=.true.
                         !print *,'equals in second word.'
                         exit
                       endif
                     enddo
                     if (equals) then
                       !Value of variable found in second word. Extract it:
                       !print *,'..found it in second word..'
                       VarValue=word(i+1:wordLength)
                       foundvariable=.true.
                       exit
                     endif
                   endif
                 endif
               endif
             endif
          endif

        endif
        deallocate(words)

       endif
    enddo
    !print *,'VarValue=',VarValue
    close(1)

end subroutine getsettingspara
      subroutine getTime(timeReturn)

      ! Return time in seconds (need to subtract the time at start of
      ! simulation to get run-time)
      !
      ! Marcel Sluiter

      implicit none

      integer nl, timetaken(6),timeReturn

      character numb*6,line*80,date*8,time*10,zone*5

      call date_and_time(date,time,zone) 

      read(date(1:4),*)timetaken(1) !years
      read(date(5:6),*)timetaken(2) !months
      read(date(7:8),*)timetaken(3) !days
      read(time(1:2),*)timetaken(4) !hours
      read(time(3:4),*)timetaken(5) !minutes
      read(time(5:6),*)timetaken(6) !seconds

      timeReturn=(((timetaken(2)*30+timetaken(3))*24+timetaken(4))*60 &
                     +timetaken(5))*60+timetaken(6)

      end subroutine getTime
subroutine initializemeam

    !--------------------------------------------------------------c
    !
    !     Initialize (initial) MEAM parameters either using a
    !     potential file provided by the user; a series of
    !     potential files in case of a continuation job; or set
    !     them to zero (later to be assigned random values) if
    !     they are to be initialized using parameters from the
    !     settings file. MEAM parameters are read both into the 
    !     work arrays: cmin, cmax, etc; and also into the p() 
    !     array.
    !
    !     Calls:         readmeamparam,variables_to_p,
    !                 p_to_variables
    !     Files read:    'startparameterfile' (filename provided
    !                 by user), potparas_best#n (n=1-10),
    !                 bestoptfuncs
    !     Returns:       cmin,cmax,meamtau,
    !                 meamrhodecay,meame0,meamrho0,meamemb3,
    !                 meamemb4,pairpotparameter,rs,rc,enconst,p()
    !
    !     Andrew Duff Jan 2015
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_atomproperties
    use m_filenames
    use m_generalinfo
    use m_optimization

    implicit none

    logical exist
    integer i,j,tmp
    character*1 tmp2
    character*20 string1
    character*80 filename,string
    

    print *
    print *,'Potential initialization'
    print *,'------------------------'

    if (.not.allocated(cmin)) then
        allocate( cmin(maxspecies,maxspecies,maxspecies), &
            cmax(maxspecies,maxspecies,maxspecies), &
            meamtau(0:lmax,maxspecies), &
            meamrhodecay(12,0:lmax,maxspecies,maxspecies), &
            meame0(maxspecies), &
            meamrho0(maxspecies), &
            meamemb3(maxspecies), &
            meamemb4(maxspecies), &
            pairpotparameter(32,maxspecies,maxspecies), &
            rs(maxspecies,maxspecies), &
            rc(maxspecies,maxspecies), &
            enconst(maxspecies)) !, &
            ! pairpotStr(splnNvals,maxspecies,maxspecies))
    endif

    cmin=0d0
    cmax=0d0
    meamtau=0d0
    meamrhodecay=0d0
    meame0=0d0
    meamrho0=0d0
    meamemb3=0d0
    meamemb4=0d0
    pairpotparameter=0d0
    rs=0d0
    rc=0d0
    enconst=0d0

    if (.not.allocated(p)) then
        allocate(p(np))
    endif

    p=0d0

    if (readpotfile.eqv..true.) then
       print *,'Reading in potential parameters from ',trim(startparameterfile)
       call readmeamparam(startparameterfile) !Read into sensibly-named variables (meamrhodecay, etc)
       call variables_to_p !Now copy these variables into the p() array (necc. for CG optimizer)
    else
       if (contjob.eqv..true.) then
          ! Continuation job: Read in potentials from files
          allocate(p_saved(np,noptfuncstore),bestoptfuncs(noptfuncstore),timeoptfunc(noptfuncstore))
          print *,'Reading in files from potparas_best#n, with #n=1,NOPTFUNCSTORE'
          open(50,file='bestoptfuncs')
          read(50,*) string
          print *,'string=',string
          if (string(1:3).ne.'Top') then
             print *,'ERROR: You do not yet have a full set of potparas_best files, stopping.'
             stop
          endif
          do i=1,noptfuncstore
             !Setup filename for each potential in turn
             filename="potparas_best"
             if (i.lt.10) then
                 write(string1,'(I1)') i
             elseif (i.lt.100) then
                 write(string1,'(I2)') i
             else
                 print *,'ERROR: more than 100 files set to save; code needs'
                 print *,'changing, STOPPING.'
                 stop
             endif
             filename=trim(filename)//trim(string1)

             !Read potential parameters from file
             print *,'starting with ',trim(adjustl(filename))
             call readmeamparam(trim(adjustl(filename)))
             call variables_to_p
             do j=1,np
                 p_saved(j,i)=p(j)
             enddo
             !Copy corresponding optimization function into bestoptfuncs
             read(50,*) tmp,tmp2,bestoptfuncs(i)
             timeoptfunc(i)=0
          enddo
          print *,'Bestopfuncs read in:'
          do i=1,noptfuncstore
             print *,i,': ',bestoptfuncs(i)
          enddo
       else
          print *,'Potential parameters to be optimized read from settings file.'
          !Set-up dummy values for now, and then after the distance analysis,
          !assign values to the cut-off radii.
          call p_to_variables
       endif
    endif

end subroutine initializemeam
subroutine getmaxspecies

    !--------------------------------------------------------------c
    !
    !     Read in the set in the vasprun.xml files, or
    !     alternatively POSCAR files, and determine
    !     the total number of species used in them.
    !     UPDATE: POSCAR support removed (relevant lines retained
    !     however so that it can be re-instated later if needed)
    !
    !     Called by:     program MEAMfit
    !     Calls:         readposcar
    !     Returns:       maxspecies 
    !
    !     Andrew Duff 2014
    !
    !--------------------------------------------------------------c

    use m_generalinfo
    use m_atomproperties
    use m_filenames
    use m_poscar
    use m_optimization
    use m_geometry
    use m_meamparameters

    implicit none
    integer i,j,nconfig_prev
    character *80 poscarfile_prev

    z2species=0
    maxspecies=0
    do i=1,nstruct !loop through each atomic configuration
        !print *,'file=',i
        if ((i.gt.1).and.(trim(poscarfiles(i)).eq.poscarfile_prev).and. &
              (nconfig(i).eq.nconfig_prev+1)) then
           !print *,'reading from left of position (nconfig(',i,')=',nconfig(i),&
           !         ', nconfig_prev+1=',nconfig_prev+1,')'
           !Do not need to read all the way from the start of the file
           call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.true. )
        else 
           !print *,'reading from start (nconfig(',i,')=',nconfig(i), &
           !          ', nconfig_prev+1=',nconfig_prev+1,')'
           !Need to read from start of the file
           call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.false. )
        endif
        !Find all atomic species
        do j=1,nspecies
            if(z2species(z(j)).eq.0) then
                maxspecies=maxspecies+1
                z2species(z(j))=maxspecies
            endif
        enddo
        poscarfile_prev=trim(poscarfiles(i))
        nconfig_prev=nconfig(i)
    enddo

end subroutine getmaxspecies





subroutine initializestruc(filenmr)

    !--------------------------------------------------------------c
    !
    !     Read in the set of vasprun.xml or POSCAR files, process 
    !     them, and store geometrical info in the module 
    !     'm_geometry'. If filenmr=0 then all files are processed, 
    !     otherwise only vasprun.xml/POSCAR file number 'filenmr' 
    !     is treated.
    !
    !     For filenmr=0, vasprun.xml/POSCAR files are read three 
    !     times, whereas for filenmr>0 only two passes are used 
    !     (and for only a single vasprun.xml/POSCAR file). The 
    !     first pass (only for filenmr=0) determines the total 
    !     number of species across all vasprun.xml/POSCAR files, 
    !     maxspecies, and constructs the array z2species, for 
    !     which the Zth element contains the order in which the 
    !     atomic number Z first appeared in the vasprun.xml/
    !     POSCAR files (and this order is henceforth refered to 
    !     as the 'species'). E.g., for a vasprun.xml file
    !     containing two species, C (Z=12) and Fe (Z=26). C is 
    !     first to appear and so we have z2species(12)=1, and 
    !     next is Fe, so that z2species(26)=2. In this way, the 
    !     z2species array maps the atomic number z to the species 
    !     number. In this example maxspecies=2.
    !
    !     The next two passes are used for both filenmr=0 and for
    !     filenmr>0, but in the latter only a single vasprun.xml/
    !     POSCAR file is processed. In the second pass, we read in 
    !     the vasprun.xml/POSCAR files again and extend the 
    !     lattice where neccesary (for cases where force 
    !     optimisation is specified) and determine the largest
    !     values of n_inequivalentsites, etc (so that we can
    !     define array sizes). Then we can read in for the third
    !     time and fill the gn_neighbors arrays, etc.
    !
    !     Notes for programmer:
    !     When 'filenmr' is non-zero, only the neighborhood of the 
    !     specified structure number is determined. However, the 
    !     arrays are still re-allocated and given sizes spanning 
    !     the full range of structure files, with array elements
    !     elements corresponding to files .ne. filenmr set to 
    !     zero. This turned out to be the tidiest way of
    !     implementing the 'filenmr' switch.
    !
    !     Written by Andrew Duff.
    !
    !     Called by:     program MEAMfit
    !     Calls:         readposcar,extendlattice,neighborhood
    !     Returns:       maxnnatoms,maxspecies,z2species,
    !                 gn_inequivalentsites,gspecies,gn_neighbors,
    !                 gneighborlist,gxyz
    !     Files read:    poscarfiles
    !     Files written: crystalstrucfile
    !
    !--------------------------------------------------------------c

    use m_filenames
    use m_poscar
    use m_neighborlist
    use m_geometry
    use m_atomproperties
    use m_datapoints
    use m_meamparameters
    use m_optimization

    implicit none
    integer i,j,maxsites,maxneighbors,ns,nshell, &
        iaux(1),filenmr,nconfig_prev,err
    character*80 poscarfile_prev

    if (verbose) then
       print *,'Preparing to initialize structures...'
    endif

    if (filenmr.eq.0) then

        !First pass:
        z2species=0
        maxspecies=0
        do i=1,nstruct
            if ((i.gt.1).and.(trim(poscarfiles(i)).eq.poscarfile_prev).and. &
              (nconfig(i).eq.nconfig_prev+1)) then
               call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.true. )
            else
               call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.false. )
            endif
            !Find all atomic species
            do j=1,nspecies
                if(z2species(z(j)).eq.0) then
                    maxspecies=maxspecies+1
                    if (maxspecies.gt.10) then
                       print *,'ERROR: maxspecies>10. Please re-size speciestoz'
                       print *,'array and then re-run.'
                       stop
                    endif
                    z2species(z(j))=maxspecies
                    speciestoz(maxspecies)=z(j)
                endif
            enddo
            poscarfile_prev=trim(poscarfiles(i))
            nconfig_prev=nconfig(i)
        enddo

        !Second pass:
        maxsites=0
        maxneighbors=0
        maxnnatoms=0
        do i=1,nstruct
           !print *,i,'/',nstruct
           !print *,'poscarfiles=',poscarfiles(i),', prev=',poscarfile_prev
           !print *,'nconfig(i)=',nconfig(i),' nconfig_prev+1=',nconfig_prev+1
            if ((i.gt.1).and.(trim(poscarfiles(i)).eq.poscarfile_prev).and. &
              (nconfig(i).eq.nconfig_prev+1)) then
               if (optimizeforce(i).gt.0) then
                  !call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.true. )
                  !TEMPORARILY COMMENTED Out BECAUSE NATOMS AND ZZ ARE CHANGED!!!
                  call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.false.)
               else
                   call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.true.)
               endif
            else
               call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.false. )
            endif
            n_inequivalentsites=natoms
            nshell=0
            if (optimizeforce(i).gt.0) then
              !print *,'extending...'
                !Need to make atoms in adjacent cells to the central cell
                !inequivalent sites, so that when the force calculation
                !takes place, these atoms are included in the energy
                !calculation.
                call extendlattice
            endif
            !Find all neighbors
            call neighborhood(optimizeforce(i))
            maxsites=max(maxsites,n_inequivalentsites)
            iaux=maxval(n_neighbors)
            maxneighbors=max(maxneighbors,iaux(1))
            maxnnatoms=max(maxnnatoms,nnatoms)
            poscarfile_prev=trim(poscarfiles(i))
            nconfig_prev=nconfig(i)
        enddo

    elseif (filenmr.gt.0) then

        !Second pass:
        call readposcar( trim(poscarfiles(filenmr)),nconfig(filenmr),vasprun(i),.false. )
        n_inequivalentsites=natoms
        !         if (filenmr.eq.23) then
        !            print *,'natoms (for file nmr 23)=',natoms
        !         endif
        nshell=0
        if (optimizeforce(filenmr).gt.0) then
            !Need to make atoms in adjacent cells to the central cell
            !inequivalent sites, so that when the force calculation
            !takes place, these atoms are included in the energy
            !calculation.
            call extendlattice
        endif
        !Find all neighbors
        call neighborhood(optimizeforce(filenmr))
        maxsites=n_inequivalentsites
        maxneighbors=maxval(n_neighbors)
        maxnnatoms=nnatoms
        !         if (filenmr.eq.23) then
        !            print *,'natoms still=',natoms
        !         endif

    endif

    !(Re)allocate relevant arrays (note: we allocate arrays
    !of size nstruct, even if filenmr>0, I.e. if we just need
    !the neighborhood of one structure. This is to avoid having to
    !rewrite other subroutines, which assume such a form for these
    !arrays).
    if (allocated(gn_inequivalentsites)) &
        deallocate(gn_inequivalentsites)
    if (allocated(gspecies)) deallocate(gspecies)
    if (allocated(gn_neighbors)) deallocate(gn_neighbors)
    if (allocated(gneighborlist)) deallocate(gneighborlist)
    if (allocated(gxyz)) deallocate(gxyz)
    if (allocated(gn_forces)) deallocate(gn_forces)
    if (allocated(gn_C)) deallocate(gn_C)
    if (allocated(optforce)) deallocate(optforce)
    if (allocated(diststr)) deallocate(diststr)
    if (allocated(dist2str)) deallocate(dist2str)
    if (allocated(dist3str)) deallocate(dist3str)
    if (allocated(dxstr)) deallocate(dxstr)
    if (allocated(dystr)) deallocate(dystr)
    if (allocated(dzstr)) deallocate(dzstr)
    allocate( gn_inequivalentsites(nstruct), &
        gspecies(maxnnatoms,nstruct), &
        optforce(maxsites,nstruct), &
        gn_neighbors(maxsites,nstruct), &
        gneighborlist(maxneighbors,maxsites,nstruct), &
        gxyz(3,maxnnatoms,nstruct), &
        gn_forces(nstruct), &
        gn_C(nstruct) )
    if (forcefit.and.fastForce) then
      allocate( diststr(maxneighbors,maxsites,0:maxsites,0:3,nstruct), &
          dist2str(maxneighbors,maxsites,0:maxsites,0:3,nstruct), &
          dist3str(maxneighbors,maxsites,0:maxsites,0:3,nstruct), &
          dxstr(maxneighbors,maxsites,0:maxsites,0:3,nstruct), &
          dystr(maxneighbors,maxsites,0:maxsites,0:3,nstruct), &
          dzstr(maxneighbors,maxsites,0:maxsites,0:3,nstruct),stat=err)
          if (err.ne.0) then
             print *,'ERROR: failed to initialized distance arrays, setting'
             print *,'FASTFORCE=false and retrying...'
             !First make sure that _none_ of them allocated
             if (allocated(diststr)) deallocate(diststr)
             if (allocated(dist2str)) deallocate(dist2str)
             if (allocated(dist3str)) deallocate(dist3str)
             if (allocated(dxstr)) deallocate(dxstr)
             if (allocated(dystr)) deallocate(dystr)
             if (allocated(dzstr)) deallocate(dzstr)
             fastForce=.false.
             allocate( diststr(maxneighbors,maxsites,0:0,0:0,nstruct), &
                 dist2str(maxneighbors,maxsites,0:0,0:0,nstruct), &
                 dist3str(maxneighbors,maxsites,0:0,0:0,nstruct), &
                 dxstr(maxneighbors,maxsites,0:0,0:0,nstruct), &
                 dystr(maxneighbors,maxsites,0:0,0:0,nstruct), &
                 dzstr(maxneighbors,maxsites,0:0,0:0,nstruct),stat=err)
          endif
    else
      allocate( diststr(maxneighbors,maxsites,0:0,0:0,nstruct), &
          dist2str(maxneighbors,maxsites,0:0,0:0,nstruct), &
          dist3str(maxneighbors,maxsites,0:0,0:0,nstruct), &
          dxstr(maxneighbors,maxsites,0:0,0:0,nstruct), &
          dystr(maxneighbors,maxsites,0:0,0:0,nstruct), &
          dzstr(maxneighbors,maxsites,0:0,0:0,nstruct),stat=err)
    endif
    if (err.ne.0) then
       print *,'ERROR: failed to initialize distance arrays, even with'
       print *,'reduced sizes (FASTFORCE=false). Stopping.'
       stop
    endif
    if (.not.allocated(gxyz_backup)) then
        allocate(gxyz_backup(3,maxnnatoms,nstruct))
        !This should just be allocated once and then left alone
    endif

    !Third pass:
    if (filenmr.eq.0) then

        do i=1,nstruct
            if ((i.gt.1).and.(trim(poscarfiles(i)).eq.poscarfile_prev).and. &
              (nconfig(i).eq.nconfig_prev+1)) then
               if (optimizeforce(i).gt.0) then
                  !call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.true. )
                  !TEMPORARILY COMMENTED Out BECAUSE NATOMS AND ZZ ARE CHANGED!!!
                  call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.false.)
               else
                  call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.true.)
               endif
            else
               call readposcar( trim(poscarfiles(i)),nconfig(i),vasprun(i),.false. )
            endif
            !Use nat and z and nspecies to calculate and then store the
            !no. of C atoms
            gn_C(i)=0
            do j=1,nspecies
                if (z(j).eq.6) then
                    gn_C(i)=nat(j)
                endif
            enddo
            n_inequivalentsites=natoms
            nshell=0
            if (optimizeforce(i).gt.0) then
                !Need to make atoms in adjacent cells to the central cell
                !inequivalent sites, so that when the force calculation
                !takes place, these atoms are included in the energy
                !calculation.
                call extendlattice
            endif
            call neighborhood(optimizeforce(i))
            ns=n_inequivalentsites
            gn_inequivalentsites(i)=ns
            !            print *,'istr=',i,', gn_inequivalentsites=',
            !     +          gn_inequivalentsites(i)
            if (optimizeforce(i).gt.0) then
                gn_forces(i)=norigatoms
            endif
            gn_neighbors(1:ns,i)=n_neighbors(1:ns)
            do j=1,ns
                gneighborlist(1:n_neighbors(j),j,i)= &
                    neighborlist(1:n_neighbors(j),j)
            enddo
            do j=1,nnatoms
                !print *,'j=',j,' species(j)=',species(j)
                gxyz(1:3,j,i)=xyz(1:3,j)!ERROR: here the second index is
                !the nearest neighbor number
                gspecies(j,i)=species(j)
            enddo
            !stop
            poscarfile_prev=trim(poscarfiles(i))
            nconfig_prev=nconfig(i)
        enddo
    else
        call readposcar( trim(poscarfiles(filenmr)),nconfig(filenmr),vasprun(i),.false. )
        !Use nat and z and nspecies to calculate and then store the
        !no. of C atoms
        gn_C(filenmr)=0
        do j=1,nspecies
            if (z(j).eq.6) then
                gn_C(filenmr)=nat(j)
            endif
        enddo
        n_inequivalentsites=natoms
        !         if (filenmr.eq.23) then
        !            print *,'n_inequiv... still=',n_inequivalentsites
        !         endif

        nshell=0
        if (optimizeforce(filenmr).gt.0) then
            !Need to make atoms in adjacent cells to the central cell
            !inequivalent sites, so that when the force calculation takes place,
            !these atoms are included in the energy calculation.
            call extendlattice
        endif
        call neighborhood(optimizeforce(filenmr))
        ns=n_inequivalentsites
        gn_inequivalentsites(filenmr)=ns
        !         if (filenmr.eq.23) then
        !            print *,'gn_inequiv(23)=',gn_inequivalentsites(23),ns
        !         endif
        if (optimizeforce(filenmr).gt.0) then
            gn_forces(filenmr)=norigatoms
        endif
        gn_neighbors(1:ns,filenmr)=n_neighbors(1:ns)
        do j=1,ns
            gneighborlist(1:n_neighbors(j),j,filenmr)= &
                neighborlist(1:n_neighbors(j),j)
        enddo
        do j=1,nnatoms
            gxyz(1:3,j,filenmr)=xyz(1:3,j)!ERROR: here the second index
            !is the nearest neighbor number
            gspecies(j,filenmr)=species(j)
        enddo
    endif
    close(1)

    if(allocated(nat)) deallocate(nat,z,coordinates)
    if(allocated(xyz)) deallocate(xyz,species, &
        n_neighbors,neighborlist)

    print *,'Completed structure initialization'

end subroutine initializestruc
subroutine initParaLimits

    !--------------------------------------------------------------c
    !
    !     Use values of 'maxspecies' and 'lmax' to set up the    
    !     variables to be used to reference the potential parameter
    !     array, p().
    !
    !     Called by:     readsettings
    !     Calls:         -
    !     Returns:       lm1, m1, m2, m3, m4
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff 2014
    !
    !--------------------------------------------------------------c

    use m_filenames
    use m_meamparameters
    use m_atomproperties

    implicit none

    lm1=lmax+1
    m1=maxspecies
    m2=m1*maxspecies
    m3=m2*maxspecies
    m4=2*m3+(4+lm1)*m1+32*m2+12*lm1*m2

end subroutine initParaLimits
subroutine map_icoords_to_its(it1,it2,it3,iface,ishell, &
        icoord1,icoord2)

    !---------------------------------------------------------------c
    !
    !     Takes the variables 'icoord1' and 'icoord2' (which are
    !     varied as:     do icoord1=-ishell,ishell
    !                       do icoord2=-ishell,ishell
    !     in the subroutine 'calculateshells') and convert them
    !     into coordinates, (it1,it2,it3), on the face of a cube.
    !     The cube is of dimension: (2 * ishell) and is centred on
    !     the origin. The value of iface denotes which face we want
    !     the coordinates to lie on.
    !
    !     Called by:     calculateshells
    !     Calls:         -
    !     Arguments:     icoord1,icoord2,iface,ishell
    !     Returns:       it1,it2,it3
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Dec 2007
    !
    !---------------------------------------------------------------c

    implicit none

    integer iface,it1,it2,it3,icoord1,icoord2,ishell

    if (iface.eq.1) then      !Map icoord1 and icoord2 onto
        it1=icoord1            !direct coordinates. iface=1 is the
        it2=icoord2            !bottom face, 2 is the top face, 3
        it3=-ishell            !the left face, 4 the right face,
    elseif (iface.eq.2) then  !5 the back face and 6 is the
        it1=icoord1            !front face.
        it2=icoord2
        it3=ishell
    elseif (iface.eq.3) then
        it1=icoord1
        it2=-ishell
        it3=icoord2
    elseif (iface.eq.4) then
        it1=icoord1
        it2=ishell
        it3=icoord2
    elseif (iface.eq.5) then
        it1=-ishell
        it2=icoord1
        it3=icoord2
    elseif (iface.eq.6) then
        it1=ishell
        it2=icoord1
        it3=icoord2
    endif

end subroutine map_icoords_to_its
subroutine matinv(a,ndim,n,b)

    !-----------------------------------------------------------c
    !
    !     Matrix inversion by Gaussian algorithm. Matrix 'a'
    !     is supplied, and it's inverse, matrix 'b' is returned.
    !     The dimension of 'a' and 'b' is ndim (=n)
    !
    !     Called by:     readposcar
    !     Calls:         -
    !     Arguments:     a,ndim,n
    !     Returns:       b
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter, April 25 2000
    !
    !-----------------------------------------------------------c

    implicit none

    integer ndim,n,i,k,imax
    real(8), parameter:: eps=1d-14
    real(8) a(ndim,ndim),b(ndim,ndim),amax,amult,div
    real(8), allocatable:: acopy(:,:),aux(:)

    allocate(acopy(n,n),aux(n))
    b=0d0
    do i=1,n
        b(i,i)=1d0
        acopy(i,1:n)=a(i,1:n)
    enddo

    do k=1,n
        amax=0d0
        imax=0
        if (k .ne. n) then
            do i=k,n
                if (abs(acopy(i,k)) .gt. amax) then
                    amax=abs(acopy(i,k))
                    imax=i
                endif
            enddo
            if (imax.ne.k) then
                aux=acopy(imax,1:n)
                acopy(imax,1:n)=acopy(k,1:n)
                acopy(k,1:n)=aux
                aux=b(imax,1:n)
                b(imax,1:n)=b(k,1:n)
                b(k,1:n)=aux
            endif
        endif
        div=acopy(k,k)
        if(abs(div).le.eps) div=eps  !dont crash if singular
        acopy(k,1:n)=acopy(k,1:n)/div
        b(k,1:n)=b(k,1:n)/div
        do i=1,n
            if (i.ne.k) then
                amult=acopy(i,k)
                acopy(i,1:n)=acopy(i,1:n)-amult*acopy(k,1:n)
                b(i,1:n)=b(i,1:n)-amult*b(k,1:n)
            endif
        enddo
    enddo
    deallocate(acopy,aux)

end subroutine matinv
subroutine meamenergy(is,iatom,cart,meam_energy,sumpp,summeamf)

    !--------------------------------------------------------------c
    !
    !     Computes the MEAM energy for structure, is.
    !
    !     Called by:     Fnew, meamforce
    !     Calls:         readmeamparam,screening_ij,
    !                 radialdensityfunction,
    !                 radialdensityfunction_onlyiatom
    !                 electrondensity,
    !                 backgrounddensity,embeddingfunction,
    !                 pairpotential,distance
    !     Arguments:     is,gn_inequivalentsites,gspecies,
    !                 gn_neighbors,gneighborlist,screening
    !     Returns:       meam_energy
    !     Files read:    none
    !     Files written: none
    !
    !     Andrew Duff and Marcel Sluiter.
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_geometry
    use m_screening   !S_ij
    use m_filenames
    use m_electrondensity!can remove
    use m_generalinfo

    implicit none

    integer i,j,jj,isp,jsp,is,iatom,cart
    real(8) pairpot,meam_energy,SI_to_eVamg, &
        summeamf,sumpp!,time!,timer

    parameter(SI_to_eVamg=14.4375d0)

    meam_energy=0d0
    sumpp=0d0

    istr=is
    call screening_ij
    if (allocated(meam_f)) deallocate(meam_f)
    allocate(meam_f(gn_inequivalentsites(istr)), &
        meam_paire(gn_inequivalentsites(istr)))

    if (iatom.eq.0) then !This is the full meamenergy calculation, and
        !calculates the total energy for the
        !supercell, for use with those datapoints
        !which require a fitting to the energy (and
        !not the forces).
        do i=1,gn_inequivalentsites(istr)
            isp=gspecies(i,istr)
            meam_paire(i)=0d0
            do jj=1,gn_neighbors(i,istr)
                j=gneighborlist(jj,i,istr)
                jsp=gspecies(j,istr)
                rij=diststr(jj,i,0,0,istr)
                !            rij=distance(i,j)
                !            if (rij.ne.diststr(jj,i,0,0,istr)) then
                !               print *,'stored dist(',diststr(jj,i,0,0,istr),
                !     +  ') differs from calculated dist (',rij,'), stopping.'
                !               print *,'jj=',jj,', i=',i,', istr=',istr,
                !     +                 ',0,0)'
                !               stop
                !            endif
                call pairpotential(isp,jsp,rij, &
                    pairpot)
                meam_paire(i)=meam_paire(i)+pairpot* &
                    screening(jj,i)
            enddo
            meam_paire(i)=meam_paire(i)*0.5d0
            sumpp=sumpp+meam_paire(i)
            meam_energy=meam_energy+meam_paire(i)
        enddo
      !  write(*,'(A9,F25.20)') 'energy(1)',meam_energy

        call radialdensityfunction

        meam_energy=meam_energy+ &
            (gn_inequivalentsites(istr)-gn_C(istr))*enconst(1) &
            +  gn_C(istr)*enconst(2)
      !  write(*,'(A9,F25.20)') 'energy(2)',meam_energy


        !       print *,'no. spc 1:',gn_inequivalentsites(istr)-gn_C(istr)
        !       print *,'no. spc 2:',gn_C(istr)

        !        print *,'meam_energy=',meam_energy
        !        print *,'enconst(1)=',enconst(1),', enconst(2)=',enconst(2), &
        !                'enconst(3)=',enconst(3)
        !        stop

    else               
        !This is for the force calculation: only those
        !embedding energies and pair-potentials which
        !have changed due to the atom being displaced
        !need to be calculated.
        isp=gspecies(iatom,istr)
        meam_paire(iatom)=0d0

        if (fastForce) then

           do jj=1,gn_neighbors(iatom,istr)
               j=gneighborlist(jj,iatom,istr)
               jsp=gspecies(j,istr)
               rij=diststr(jj,iatom,iatom,cart,istr)
               !            rij=distance(iatom,j)
               !            if (rij.ne.diststr(jj,iatom,iatom,cart,istr)) then
               !              print *,'stored dist(',diststr(jj,iatom,iatom,cart,istr),
               !     +  ') differs from calculated dist (',rij,'), stopping.'
               !              print *,'jj=',jj,', iatom=',iatom,', istr=',istr,
               !     +                ', iatom=',iatom,', cart=',cart
               !              stop
               !            endif
               call pairpotential(isp,jsp,rij, &
                   pairpot)

               meam_paire(iatom)=meam_paire(iatom)+pairpot* &
                   screening(jj,iatom)
           enddo

        else

           do jj=1,gn_neighbors(iatom,istr)
               j=gneighborlist(jj,iatom,istr)
               jsp=gspecies(j,istr)
               rij=distance(iatom,j)
               call pairpotential(isp,jsp,rij, &
                   pairpot)

               meam_paire(iatom)=meam_paire(iatom)+pairpot* &
                   screening(jj,iatom)
           enddo

        endif

        sumpp=meam_paire(iatom)
        meam_energy=meam_energy+meam_paire(iatom)

        call radialdensityfunction_onlyiatom(iatom,cart)

    endif

    call electrondensity(iatom,cart)
    call backgrounddensity
    call embeddingfunction

    ! write(*,'(A9,F25.20)') 'energy(2B)',meam_energy
    summeamf=0d0
    do i=1,gn_inequivalentsites(istr)
        !meam_energy=meam_energy+meam_f(i)
        summeamf=summeamf+meam_f(i)
    enddo
    meam_energy = meam_energy + summeamf

   !if (istr.eq.1) then
   !  write(*,'(A6,F25.20)') 'energy',meam_energy
   !  write(*,'(A6,F25.20,A11,F30.20)') 'sumpp=',sumpp,', summeamf=',summeamf
   !  write(*,'(A6,F30.20,A11,F30.20)') 'enZrpt',(gn_inequivalentsites(istr)-gn_C(istr))*enconst(1),  &
   !                                    'enC        ',gn_C(istr)*enconst(2)
   !   print *,'tot_enconst=',(gn_inequivalentsites(istr)-gn_C(istr))*enconst(1) &
   !        +  gn_C(istr)*enconst(2)
   !    print *,'TOTAL=',meam_energy
   !endif

    deallocate(meam_t,meam_f,meam_paire)

end subroutine meamenergy
 

program MEAMfit 

    !----------------------------------------------------------------------c
    !
    !     "MEAMfit" version 1.01, written by Andrew I. Duff and Marcel H. F.
    !     Sluiter, Copyright, 2015.
    !
    !     This program optimizes MEAM parameters so that the MEAM energies
    !     of a particular set of crystal structures and the forces acting on
    !     the atoms in these structures are equal to the same quantities as
    !     calculated using the VASP DFT program.
    !
    !     A full description of the program as well as citation details
    !     is given in the MEAMfit User Guide.
    !
    !     Abbreviations used in commenting: optfunc=optimization function;
    !     CH=conjugate gradient; GA=genetic algorithm
    !
    !     When using this program please cite the following article:
    !     ...
    !
    !     LEGAL NOTE:
    !
    !     The authors make no warranty, express or implied, with respect
    !     to this software, its quality, performance, merchantability,
    !     or fitness for a particular purpose. In no event shall anyone
    !     involved in the creation of this software, be liable for direct,
    !     indirect, special, or consequential damages. The software is
    !     made available "as is", and the user assumes the entire risk of
    !     its performance and suitability to the user's purpose.
    !
    !     The code is free to use and distribute, however it must be
    !     provided in its original, non-modified form. Any requests for 
    !     modification should be directed to Andrew I. Duff (see MEAMfit 
    !     manual for latest contact details).
    !
    !----------------------------------------------------------------------c

    use m_optimization
    use m_filenames
    use m_geometry
    use m_datapoints
    use m_generalinfo
    use m_meamparameters

    implicit none
    logical maxOK
    integer i,nfreep,iseed,nattempts
    real(8) tmp

    call getTime(tstart)

    print *
    print *,'-------------- MEAMfit (version 1.01) ---------------'
    print *,'By Andrew I. Duff and Marcel H. F. Sluiter, 2006-2015'
    print *,'-----------------------------------------------------'

    call testargs !Check for command-line arguments
    call setupfilenames  !Define all filenames
    call checkInputFilesExist !Check if 'fitdbse' and 'settings' files exist
    call getmaxspecies   !Determine 'maxspecies' from the vasprun files (the total number of species
                  !to be used). This needs doing in advance of readsettings for setting up arrays.
    if (readpotfile.eqv..true.) then
        call extractlmax
    endif
    if (settingsfileexist.eqv..false.) then
       call setupdefaultlmax
    endif

    if (settingsfileexist.eqv..true.) then
       print *
       print *,'General initialization'
       print *,'----------------------'
    endif
 
    call readsettings(nfreep,iseed) !read settings from settings file
    if (settingsfileexist.eqv..true.) then
       call setuprandomseed(iseed)
    endif
    !---- Initial meam-parameter set-up ----
    call initializemeam       !Initialize initial meam parameters
    if ((readParasOnebyone.eqv..false.).and.(readpotfile)) then
       !If a potential file is supplied by user, amend freep so that non-zero parameters
       !in this file are optimized (unless PARASTOOPT is specified in the settings file)
       call amendFreep
    endif
    call getRmax !Find the maximum interatomic separation that needs to be considered.

  ! !Prepare splines
  ! narrpairpotXX=1000!1000
  ! narrthiXX=1000
  ! allocate(r_pairpotXXarr(narrpairpotXX), &
  !     pairpotXXarr(narrpairpotXX), &
  !     secderpairpotXX(narrpairpotXX), &
  !     r_thiXXarr(narrthiXX), &
  !     thiXXarr(narrthiXX), &
  !     secderthiXX(narrthiXX))
  ! do i=1,narrpairpotXX
  !     r_pairpotXXarr(i) = &
  !         ( (dble(i)/dble(narrpairpotXX)) )**2 * &
  !         p_rmax
  ! enddo
  ! do i=1,narrthiXX
  !     r_thiXXarr(i) = &
  !         ( (dble(i)/dble(narrthiXX)) )**2 * &
  !         p_rmax
  ! enddo
    !--------------------------------------

    if (optimizestruc) then
        !Optimize positons of ions, keeping the MEAM parameters fixed.
        !Note: currently unsupported - to be fully instated later
        call optimizestructure(.false.) !Index: optimize all structures?
        stop
    endif

    !---- Set-up structures and fitting database ----
    call initializestruc(0) !Read in atomic configurations from vasprun files
    call backupxyz !Store original positions so we can reset them after each optimization
    !cycle.
    call readdatapoints  !Read energies and forces from vasprun files
    call calcSds       !Determine s.d's of energies and forces to rescale optfunc
    call setupnntables !Set up atom separation tables, also for forces (must
                       !come after 'readdatapoints' as this determines which
                       !forces are to be optimized)

    call findsmallestsepn !Determine smallest and largest interatomic separations
    if ((startparas.eq.1).and.(readpotfile.eqv..false.)) then
       print *,'sepnHistogram file produced. No potential filename supplied and'
       print *,'-noopt/-no flag at command-line, STOPPING.'
       stop
    endif

    !------------------------------------------------

    if ((settingsfileexist.eqv..false.).or. &
        (writepotfile.eqv..true.)) then
       !Create template file for 'potparas_best' if requested
       call createMeamParasTemplate
    endif

    if (startparas.eq.2) then
       !Parameters to be randomly initialized: Read in bounds for initialization
       if (readparasfromsettings) then
          call readBoundsForRandParas
          print *,"Read in limits on random parameters from settings file"
       else
          call defaultBoundsForRandParas
          print *,"Limits on random parameters set to default values"
       endif
       print *
       print *,"Limits used to initialize potential parameters:"
       call displayBoundsForRandParas
    endif

    if (printoptparas.eqv..true.) then
       !If requested, write out sample PARASTOOPT and stop
       call writeoptparasSub
       print *
       print *,'Stopping.'
       stop
    endif
 
    if (nsteps.gt.1) then
       call checkParas !Check bounds for random parameter initialization as well
                       !as the parameters themselves if supplied by user
       print *
       print *,'Beginning optimization'
       print *,'----------------------'
       print *
    endif

    if (contjob.eqv..false.) then
       n_optfunc=1
       allocate(bestoptfuncs(noptfuncstore),timeoptfunc(noptfuncstore))
    else
       n_optfunc=noptfuncstore+1
    endif

    !Main optimization loop: randomly seed MEAM parameters; optimize using
    !conjugate-gradient; repeat until optimization function small enough.
    lowestoptfuncrand=100000d0
    do

        if (startparas.eq.2) then
           call generaterandomstartparas !Randomly initialize starting potential parameters
           if ((n_optfunc.gt.noptfuncstore).and.(genAlgo)) then
              !Once we have full set of parameters, start genetic algorithm:
              print *,'Mix potential from existing potentials; with mutations.'
              call mixPotential
           endif
           call variables_to_p !Move parameters to the p() array (for CG optimization)
        endif
 
        call optimizeparameters !Main call: Optimize the parameters

        if (startparas.eq.1) then
            print *,'PRAXIS failed to find a minimum, and'
            print *,'startparas set to 1 therefore'
            print *,'stopping.'
            print *,'Waiting for cycle restart...'
            stop
        endif

        if ((opt_failed.eqv..false.).and.(readpotfile.eqv..true.)) then
           !Some of the (sensible-named) MEAM parameters read in from file will have been
           !changed in the optimization. Reset them here.
           call readmeamparam(startparameterfile)
        endif

    enddo

end program MEAMfit
subroutine meamforce(i,iatom,force)

    !------------------------------------------------------------c
    !
    !     Calculate the force 'force' acting on atom 'iatom' of
    !     structure 'i'.
    !
    !     The method is to evaluate the meamenergy when atom
    !     iatom is in the positions:
    !     (gxyz(1,iatom,i),gxyz(2,iatom,i),gxyz(3,iatom,i)) and
    !     rnew(1:3)=(gxyz(1,iatom,i)+delta,gxyz(2,iatom,i),
    !     gxyz(3,iatom,i)). The first energy is subtracted form
    !     the second and the overall quantity is divided by delta
    !     to give the x-component of the force. The y and z
    !     components are obtained by moving the '+delta' from the
    !     x component to the y and z components respectively.
    !
    !     Called by:     optimizationfunction
    !     Calls:         meamenergy
    !     Arguments:     i,iatom,gxyz
    !     Returns:       force
    !     Files read:    none
    !     Files written: none
    !
    !     Andrew Duff Dec 2007
    !
    !------------------------------------------------------------c

    use m_generalinfo
    use m_geometry
    use m_electrondensity

    implicit none

    integer i,iatom
    real(8) forcenew(3)
    real(8) force(3),meam_energy,meam_energy2,sumpp,summeamf
    real(8), parameter:: dist=1d-8 !The distance moved in the x,y and
    !z directions to generate the point

    !The radial densities between atoms i and j (where i.ne.iatom and
    !j.ne.iatom) need only be calculated once. For simplicity,
    !calculate the radial densities for all atoms, and then later
    !when we move atom iatom, just recalculate the radial densities
    !involving atom iatom
    istr=i
    call radialdensityfunction

    !First evaluate the meam energy for structure 'i' with atom
    !'iatom' in it's true position
    call meamenergy(i,iatom,0,meam_energy,sumpp,summeamf)
    !Now displace the coordinates of atom 'iatom' in the x-direction
    !and repeat.
    gxyz(1,iatom,i)=gxyz(1,iatom,i)+dist
    call meamenergy(i,iatom,1,meam_energy2,sumpp,summeamf)

    !We can now evaluate the x-component of the force on this atom.
    force(1)=-(meam_energy2-meam_energy)/dist
    gxyz(1,iatom,i)=gxyz(1,iatom,i)-dist !Return the atom to it's
    !correct position
    call rmforcenoise(meam_energy,meam_energy2,dist,force(1), &
       forcenew(1)) !Place zeros in the digits for the force where rounding-error occurs

    !Repeat for y component
    gxyz(2,iatom,i)=gxyz(2,iatom,i)+dist
    call meamenergy(i,iatom,2,meam_energy2,sumpp,summeamf)

    force(2)=-(meam_energy2-meam_energy)/dist
    gxyz(2,iatom,i)=gxyz(2,iatom,i)-dist
    call rmforcenoise(meam_energy,meam_energy2,dist,force(2), &
        forcenew(2))

    !Repeat for z component
    gxyz(3,iatom,i)=gxyz(3,iatom,i)+dist
    call meamenergy(i,iatom,3,meam_energy2,sumpp,summeamf)
    force(3)=-(meam_energy2-meam_energy)/dist
    gxyz(3,iatom,i)=gxyz(3,iatom,i)-dist
    call rmforcenoise(meam_energy,meam_energy2,dist,force(3), &
        forcenew(3))

    !     write(*,'(A16,F30.20,F30.20,F30.20)') 'original force:',force(1), &
    !                                           force(2),force(3)
    ! Use the following line if using rmforcenoise:
     force(1:3)=forcenew(1:3)
    !   write(*,*) ' force(i=',i,',iatom=',iatom,')=     ',force(1),force(2),force(3)


    !     Old code using three points per force component - needs updating if we
    !     want to use.
    !
    !     !New code computes both the force using a displacement in the
    !     !positive and negative directions of x, y or z, then takes the
    !     !average. IMPORTANT: if we want to start using the method below
    !     !again, need to update the meamenergy calls so they include
    !     !cartesian component, and generalise code so a minus cart comp
    !     !corresponds to a displacement in neg. dir.
    !
    !     !First evaluate the meam energy for structure 'i' with atom
    !     !'iatom' in it's true position
    !     call meamenergy(i,iatom,meam_energy,sumpp,summeamf)
    !     !Now displace the coordinates of atom 'iatom' in the x-direction
    !     !and repeat (for both +ve and -ve x directions)
    !     gxyz(1,iatom,i)=gxyz(1,iatom,i)+dist
    !     call meamenergy(i,iatom,meam_energy2,sumpp,summeamf)
    !     gxyz(1,iatom,i)=gxyz(1,iatom,i)-2d0*dist
    !     call meamenergy(i,iatom,meam_energy3,sumpp,summeamf)
    !
    !     !We can now evaluate the x-component of the force on this atom.
    !     force(1)=-0.5d0*(meam_energy2-meam_energy3)/dist
    !     gxyz(1,iatom,i)=gxyz(1,iatom,i)+dist !Return the atom to it's
    !                                          !correct position
    !     !Repeat for y component
    !     gxyz(2,iatom,i)=gxyz(2,iatom,i)+dist
    !     call meamenergy(i,iatom,meam_energy2,sumpp,summeamf)
    !     gxyz(2,iatom,i)=gxyz(2,iatom,i)-2d0*dist
    !     call meamenergy(i,iatom,meam_energy3,sumpp,summeamf)
    !     force(2)=-0.5d0*(meam_energy2-meam_energy3)/dist
    !     gxyz(2,iatom,i)=gxyz(2,iatom,i)+dist
    !
    !     !Repeat for z component
    !     gxyz(3,iatom,i)=gxyz(3,iatom,i)+dist
    !     call meamenergy(i,iatom,meam_energy2,sumpp,summeamf)
    !     gxyz(3,iatom,i)=gxyz(3,iatom,i)-2d0*dist
    !     call meamenergy(i,iatom,meam_energy3,sumpp,summeamf)
    !     force(3)=-0.5d0*(meam_energy2-meam_energy3)/dist
    !     gxyz(3,iatom,i)=gxyz(3,iatom,i)+dist

    !      print *,'force:',force(1:3)
    !      print *,'after removal of numerical noise, force:',forcenew(1:3)
    !      stop

end subroutine meamforce

       subroutine mixPotential

       ! Randomly selects two potentials out of those in the high-score
       ! board and mixs them together to create a new trial potential.
       !
       ! Andy Duff, 2014

       use m_optimization
       use m_meamparameters
       use m_atomproperties

       implicit none

       integer pot1,pot2,i,j,k,l,m
       real(8) tmp

       ! Select two potentials at random
       call RANDOM_NUMBER(tmp)
       pot1=tmp*dble(noptfuncstore)+1
       do
          call RANDOM_NUMBER(tmp)
          pot2=tmp*dble(noptfuncstore)+1
     !    print *,pot1,pot2
          if (pot1.ne.pot2) exit
       enddo
     !print *,'Preparing to mix pot ',pot1,' and pot ',pot2
     !print *
     !print *,'pot1:'
     !print *,p_saved(1:np,pot1)
     !print *
     !print *,'pot2:'
     !print *,p_saved(1:np,pot2)
     !print *
     !print *,'new pot:'

!        do j=min(n_optfunc-1,noptfuncstore-1),i,-1
!            bestoptfuncs(j+1)=bestoptfuncs(j)
!            do k=1,np
!                p_saved(k,j+1)=p_saved(k,j)
!            enddo
!        enddo

       j=1
       k=1
       l=1
       do i=1,m3
           call RANDOM_NUMBER(tmp)
           if (tmp.gt.probRand) then !Otherwise leave cmin, cmax to
                            !their previous, randomly generated, values.
              call RANDOM_NUMBER(tmp)
              if (tmp.lt.probPot1) then
                 cmin(j,k,l)=p_saved(i,pot1)
              else
                 cmin(j,k,l)=p_saved(i,pot2)
              endif
           endif
           !cmin(j,k,l)=p(i)
           if (j.lt.maxspecies) then
               j=j+1
           elseif (k.lt.maxspecies) then
               j=1
               k=k+1
           elseif (l.lt.maxspecies) then
               j=1
               k=1
               l=l+1
           endif
       enddo
       ! print *,'cmin=',cmin

       j=1
       k=1
       l=1
       do i=m3+1,2*m3
           call RANDOM_NUMBER(tmp)
           if (tmp.gt.probRand) then
              call RANDOM_NUMBER(tmp)
              if (tmp.lt.probPot1) then
                 cmax(j,k,l)=p_saved(i,pot1)
              else
                 cmax(j,k,l)=p_saved(i,pot2)
              endif
           endif
           !cmax(j,k,l)=p(i)
           if (j.lt.maxspecies) then
               j=j+1
           elseif (k.lt.maxspecies) then
               j=1
               k=k+1
           elseif (l.lt.maxspecies) then
               j=1
               k=1
               l=l+1
           endif
       enddo
       ! print *,'cmax=',cmax

       j=1
       k=0
       do i=2*m3+1,2*m3+lm1*m1
          call RANDOM_NUMBER(tmp)
          if (tmp.gt.probRand) then
             call RANDOM_NUMBER(tmp)
             if (tmp.lt.probPot1) then
                meamtau(k,j)=p_saved(i,pot1)
             else
                meamtau(k,j)=p_saved(i,pot2)
             endif
          endif
          !meamtau(k,j)=p(i)
          if (j.lt.maxspecies) then
              j=j+1
          elseif (k.lt.lmax) then
              j=1
              k=k+1
          endif
       enddo
       ! print *,'meamtau=',meamtau

       j=1
       k=0
       l=1
       m=1
       do i=2*m3+lm1*m1+1,2*m3+lm1*m1+12*lm1*m2
           if (j.eq.1) then
              call RANDOM_NUMBER(tmp)
              if (tmp.gt.probRand) then
                 call RANDOM_NUMBER(tmp)
              else
                 tmp=0d0 !Use this to signal mutation here
                 print *,'rand gen a density'
              endif
           endif
           if (tmp.ne.0d0) then
              if (tmp.lt.probPot1) then
                 meamrhodecay(j,k,l,m)=p_saved(i,pot1)
              else
                 meamrhodecay(j,k,l,m)=p_saved(i,pot2)
              endif
           endif
           !meamrhodecay(j,k,l,m)=p(i)
           if (j.lt.12) then
               j=j+1
           elseif (k.lt.lmax) then
               j=1
               k=k+1
           elseif (m.lt.maxspecies) then
               j=1
               k=0
               m=m+1
           elseif (l.lt.maxspecies) then
               j=1
               k=0
               l=l+1
               m=1
           endif
       enddo
       ! print *,'meamrhodecay=',meamrhodecay

       call RANDOM_NUMBER(tmp)
       if (tmp.gt.probRand) then
          call RANDOM_NUMBER(tmp)
          if (tmp.lt.probPot1) then
             meame0(1:maxspecies)=p_saved(2*m3+lm1*m1+12*lm1*m2+1: &
                 2*m3+m1+lm1*m1+12*lm1*m2,pot1)
          else
             meame0(1:maxspecies)=p_saved(2*m3+lm1*m1+12*lm1*m2+1: &
                 2*m3+m1+lm1*m1+12*lm1*m2,pot2)
          endif
       else
       !   print *,'rand gen an emb para'
       endif
       ! print *,'meame0=',meame0

       call RANDOM_NUMBER(tmp)
       if (tmp.gt.probRand) then
          call RANDOM_NUMBER(tmp)
          if (tmp.lt.probPot1) then
            meamrho0(1:maxspecies)=p_saved(2*m3+m1+lm1*m1+12*lm1*m2+1: &
             2*m3+2*m1+lm1*m1+12*lm1*m2,pot1)
          else
            meamrho0(1:maxspecies)=p_saved(2*m3+m1+lm1*m1+12*lm1*m2+1: &
             2*m3+2*m1+lm1*m1+12*lm1*m2,pot2)
          endif
       endif
       ! print *,'meamrho0=',meamrho0

       call RANDOM_NUMBER(tmp)
       if (tmp.gt.probRand) then
         call RANDOM_NUMBER(tmp)
         if (tmp.lt.probPot1) then
          meamemb3(1:maxspecies)=p_saved(2*m3+2*m1+lm1*m1+12*lm1*m2+1: &
              2*m3+3*m1+lm1*m1+12*lm1*m2,pot1)
         else
          meamemb3(1:maxspecies)=p_saved(2*m3+2*m1+lm1*m1+12*lm1*m2+1: &
              2*m3+3*m1+lm1*m1+12*lm1*m2,pot2)
         endif
       endif
       ! print *,'meamemb3=',meamemb3

       call RANDOM_NUMBER(tmp)
       if (tmp.lt.probPot1) then
          meamemb4(1:maxspecies)=p_saved(2*m3+3*m1+lm1*m1+12*lm1*m2+1: &
              2*m3+4*m1+lm1*m1+12*lm1*m2,pot1)
       else
          meamemb4(1:maxspecies)=p_saved(2*m3+3*m1+lm1*m1+12*lm1*m2+1: &
              2*m3+4*m1+lm1*m1+12*lm1*m2,pot2)
       endif
       !meamemb4(1:maxspecies)=p(2*m3+3*m1+lm1*m1+12*lm1*m2+1: &
       !    2*m3+4*m1+lm1*m1+12*lm1*m2)
       ! print *,'meamemb4=',meamemb4

       j=1
       if (maxspecies.eq.1) then
           call RANDOM_NUMBER(tmp)
           if (tmp.gt.probRand) then
              call RANDOM_NUMBER(tmp)
           else
              tmp=0d0 !Use this to signal mutation here
            !  print *,'rand genning a pairpot para'
           endif
           do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,2*m3+(4+lm1)*m1+12*lm1*m2+32
               if (tmp.ne.0d0) then
                  if (tmp.lt.probPot1) then
                     pairpotparameter(j,1,1)=p_saved(i,pot1)
                  else
                     pairpotparameter(j,1,1)=p_saved(i,pot2)
                  endif
               endif
               !pairpotparameter(j,1,1)=p(i)
               j=j+1
           enddo
       else
           k=1
           l=1
           do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,m4
               if (j.eq.1) then
                  call RANDOM_NUMBER(tmp)
                  if (tmp.gt.probRand) then
                     call RANDOM_NUMBER(tmp)
                  else
                     tmp=0d0 !Use this to signal mutation here
                  !   print *,'rand genning a pairpot para'
                  endif
               endif
               if (tmp.ne.0d0) then
                  if (tmp.lt.probPot1) then
                     pairpotparameter(j,k,l)=p_saved(i,pot1)
                  else
                     pairpotparameter(j,k,l)=p_saved(i,pot2)
                  endif
               endif
               !pairpotparameter(j,k,l)=p(i)
               if (j.lt.32) then
                   j=j+1
               elseif (l.lt.maxspecies) then
                   j=1
                   l=l+1
               elseif (k.lt.maxspecies) then
                   j=1
                   l=1
                   k=k+1
               endif
           enddo
       endif
       ! print *,'pairpotparameter=',pairpotparameter

       j=1
       k=1
       do i=m4+1,m4+m2
           call RANDOM_NUMBER(tmp)
           if (tmp.gt.probRand) then
              call RANDOM_NUMBER(tmp)
              if (tmp.lt.probPot1) then
                 rs(j,k)=p_saved(i,pot1)
              else
                 rs(j,k)=p_saved(i,pot2)
              endif
           endif
           !rs(j,k)=p(i)
           ! print *,'rs(',j,',',k,')=',rs(j,k)
           if (j.lt.maxspecies) then
               j=j+1
           elseif (k.lt.maxspecies) then
               j=1
               k=k+1
           endif
       enddo

       j=1
       k=1
       do i=m4+m2+1,m4+2*m2
           call RANDOM_NUMBER(tmp)
           if (tmp.gt.probRand) then
              call RANDOM_NUMBER(tmp)
              if (tmp.lt.probPot1) then
                 rc(j,k)=p_saved(i,pot1)
              else
                 rc(j,k)=p_saved(i,pot2)
              endif
           endif
           !rc(j,k)=p(i)
           ! print *,'rc(',j,',',k,')=',rc(j,k)
           if (j.lt.maxspecies) then
               j=j+1
           elseif (k.lt.maxspecies) then
               j=1
               k=k+1
           endif
       enddo

       j=1
       do i=m4+2*m2+1,m4+2*m2+m1
           call RANDOM_NUMBER(tmp)
           if (tmp.gt.probRand) then
              call RANDOM_NUMBER(tmp)
              if (tmp.lt.probPot1) then
                 enconst(j)=p_saved(i,pot1)
              else
                 enconst(j)=p_saved(i,pot2)
              endif
           else
           !   print *,'rand gening an enconst (=',enconst,')'
           endif
           !enconst(j)=p(i)
           ! print *,'enconst(',j,')=',enconst(j)
           j=j+1
       enddo

       !stop

       end subroutine mixPotential 

subroutine neighborhood(forcecalc)

    !----------------------------------------------------------------------c
    !
    !     Determine which atoms in the system are nearest neighbors (nn's)
    !     to at least one of the central atoms. Label each of these atoms
    !     from 1 to 'nnatoms' (the total number of such atoms), giving the
    !     central atoms the labels 1, ..., 'natoms' (Note: atoms are nn's
    !     if they are within 'p_rmax' of one another). Record the positions
    !     and species of these atoms in the arrays: 'xyz(1:3,1:nnatoms)'
    !     and 'species(1:nnatoms)'.
    !
    !     For each central atom, i, record which atoms are it's nn's in the
    !     array: 'neighborlist(1:n_neighbors(i),i)' (where 'n_neighbors(i)'
    !     is number of nearest neighbors of atom i). The elements of this
    !     array are just the labels defined in the previous paragraph.
    !
    !     For force-calculations: Note that the lattice has already been
    !     extended by the subroutine 'extendlattice' so that the
    !     atoms within sepn <= p_rmax of atoms in the original unit-cell
    !     are defined as inequivalent atoms (central atoms).
    !
    !     Called by:     initializestructuresandmeamparas
    !     Calls:         calculateshells
    !     Arguments:     natoms,nspecies,nat,crystalstrucfile,z,z2species,
    !                 coordinates,p_rmax,tr
    !     Returns:       species,xyz,nnatoms,n_neighbors,species,
    !                 neighborlist,n_inequivalentsites
    !     Files read:    filein
    !     Files written: crystalstrucfile
    !
    !     Marcel Sluiter, Feb 8 2006
    !     continued by Andy Duff, Dec, 2007
    !
    !----------------------------------------------------------------------c

    use m_filenames
    use m_generalinfo
    use m_atomproperties
    use m_poscar
    use m_neighborlist
    use m_datapoints

    implicit none

    logical plot_crystal_struc !Set to true to plot all nn's to a file:
                               !crystal_struc.xyz
    integer i,it1,it2,it3,j,nshell, &
        nnatoms_tot,atomnumber,forcecalc
    real(8) rmax2,xyztry1(3),xyztry2(3),xyztry(3),dxyz(3),dd
    logical inlist

    plot_crystal_struc=.false.

    if(allocated(xyz)) deallocate(xyz,species, &
        n_neighbors,neighborlist)

    n_inequivalentsites=natoms
    !Determine how many shells of unit cells should be included in
    !order that all nn's of the central atoms are included in the
    !calculation.

    if (forcecalc.gt.0) then
        call calculateshells_forces(nshell,nnatoms_tot,p_rmax)
    else
        call calculateshells(nshell,nnatoms_tot,p_rmax)
    endif
    allocate(species(nnatoms_tot),xyz(3,nnatoms_tot), &
        n_neighbors(natoms),neighborlist(nnatoms_tot,natoms))

    if (plot_crystal_struc) then
        !Plot the crystal structure
        !Note: old code which needs updating. Currently not supported.
        print *,'ERROR: plot_crystal_struc currently unsupported, stopping.'
        stop
        open(unit=50,file=trim(crystalstrucfile),status='unknown')
        write(50,*) nnatoms_tot
        write(50,*) 'Cu'
    endif

    n_inequivalentsites = natoms
    do i=1,n_inequivalentsites
        !print *,'atom ',i,' zz(i)=',zz(i),' z2species(zz(i))=',z2species(zz(i))
        species(i)=z2species(zz(i))
        xyz(1:3,i)=coordinates(1:3,i)
        if (plot_crystal_struc) then
            if ((forcecalc.gt.0).and.(i.le.natoms_old)) then
                if (species(i).eq.1) then
                    write(50,*) 'Ce',xyz(1:3,i)
                elseif (species(i).eq.2) then
                    write(50,*) 'Ta',xyz(1:3,i)
                elseif (species(i).eq.3) then
                    write(50,*) 'Mo',xyz(1:3,i)
                else
                    print *,'(1)',species(i),i
                    stop
                endif
            else
                if (species(i).eq.1) then
                    write(50,*) 'H',xyz(1:3,i)
                elseif (species(i).eq.2) then
                    write(50,*) 'Ta',xyz(1:3,i)
                elseif (species(i).eq.3) then
                    write(50,*) 'P',xyz(1:3,i)
                else
                    print *,'(1)',species(i),i
                    stop
                endif
            endif
        endif
    enddo
    !stop
    if (forcecalc.eq.0) then

        rmax2=p_rmax**2

        !Determine which of the atoms in the trial unit cells are
        !within p_rmax of the atoms in the central unit cell.
        nnatoms=n_inequivalentsites
        n_neighbors=0
        do i=1,natoms
            do it1=-nshell,nshell
                xyztry1=xyz(1:3,i)+tr(1:3,1)*real(it1)
                do it2=-nshell,nshell
                    xyztry2=xyztry1+tr(1:3,2)*real(it2)
                    do it3=-nshell,nshell
                        xyztry=xyztry2+tr(1:3,3)*real(it3)
                        inlist=(it1.eq.0.and.it2.eq.0.and.it3.eq.0)
                        atomnumber=i
                        do j=1,n_inequivalentsites
                            dxyz=xyztry-xyz(1:3,j)
                            dd=dxyz(1)**2+dxyz(2)**2+dxyz(3)**2
                            if(dd.le.rmax2) then
                                if(.not.inlist) then
                                    nnatoms=nnatoms+1 !Add to the list
                                    atomnumber=nnatoms
                                    xyz(1:3,atomnumber)=xyztry(1:3)
                                    species(atomnumber)=species(i)
                                    inlist=.true.
                                    if (plot_crystal_struc) then
                                        if (species(atomnumber).eq.1) then
                                            write(50,*) 'Mg',xyztry
                                        elseif (species(atomnumber).eq.2) then
                                            write(50,*) 'Ca',xyztry
                                        elseif (species(atomnumber).eq.3) then
                                            write(50,*) 'P',xyztry
                                        else
                                            print *,'!',species(atomnumber)
                                            stop
                                        endif
                                    endif
                                endif
                                if(atomnumber.ne.j) then
                                    n_neighbors(j)=n_neighbors(j)+1 !Update
                                    !neighborlist
                                    neighborlist(n_neighbors(j),j)=atomnumber
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

    else

        rmax2=p_rmax**2
        !Determine which of the atoms in the trial unit cells are within
        !p_rmax of the inequivalent atoms.
        nnatoms=n_inequivalentsites
        n_neighbors=0
        do i=1,natoms_old
            do it1=-nshell,nshell
                xyztry1=coords_old(1:3,i)+tr(1:3,1)*real(it1)
                do it2=-nshell,nshell
                    xyztry2=xyztry1+tr(1:3,2)*real(it2)
                    do it3=-nshell,nshell
                        xyztry=xyztry2+tr(1:3,3)*real(it3)
                        call checkatom(xyztry,inlist,atomnumber)
                        do j=1,n_inequivalentsites
                            dxyz=xyztry-xyz(1:3,j)
                            dd=dxyz(1)**2+dxyz(2)**2+dxyz(3)**2
                            if(dd.le.rmax2) then
                                if(.not.inlist) then
                                    nnatoms=nnatoms+1 !Add to the list
                                    atomnumber=nnatoms
                                    xyz(1:3,atomnumber)=xyztry(1:3)
                                    species(atomnumber)=species(i)
                                    inlist=.true.
                                    if (plot_crystal_struc) then
                                        if (species(atomnumber).eq.1) then
                                            write(50,*) 'Mg',xyztry
                                        elseif (species(atomnumber).eq.2) then
                                            write(50,*) 'Ta',xyztry
                                        elseif (species(atomnumber).eq.3) then
                                            write(50,*) 'P',xyztry
                                        else
                                            print *,'!',species(atomnumber)
                                            stop
                                        endif
                                    endif
                                endif
                                if(atomnumber.ne.j) then
                                    n_neighbors(j)=n_neighbors(j)+1 !Update
                                    !neighborlist
                                    neighborlist(n_neighbors(j),j)=atomnumber
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

    endif

    if (plot_crystal_struc) then
        !Plot the crystal structure
        close(unit=50)
        print *,'Finished plotting crystal structure.'
        print *,'Stopping.'
        stop
    endif

end subroutine neighborhood
subroutine optimizeparameters

    !----------------------------------------------------------------------c
    !
    !      Optimizes parameters defined in p().
    !
    !      Called by:     program MEAMfit
    !      Calls:         optimizationfunction,random_number,
    !                 writemeamp,displayparas
    !      Returns:       bestfitdata
    !      Files read:    -
    !      Files written: -
    !
    !      Andrew Duff 2010-2015
    !
    !----------------------------------------------------------------------c

    use m_datapoints
    use m_filenames
    use m_optimization
    use m_generalinfo
    use m_atomproperties
    use m_plotfiles
    use m_geometry

    use m_meamparameters

    implicit none

    logical largest,signflipped(np)
    integer hrs
    character*80 filename
    character*20 string1,string2,string3
    integer, parameter:: pmaxtry=1
    integer, parameter:: nPotsSave=3
    real(8), parameter:: tolerance=1d-9,maxchange=1d0
    integer i,j,k,ip,iter,nfreep,n_prx,tmpint,lv,nf
    real(8) func,func0,bestfunc,aux1,aux2,p_old,y2(21),tmp,tmp2, &
        praxis,fmin,cutoff(np),func_prev,func_diff, &
        prx_thresh
    real(8), allocatable:: popt(:),d_smsno(:),urp(:,:)
    !Variables for CG optimizer:
    integer uip(1)
    integer, parameter :: liv=60
    integer iv(liv)
    real(8), allocatable:: vv(:)

    external Fnew

    !Initialize P
    if (.not.allocated(positivep)) then
        allocate(scaledstep(np),failedtries(np),positivep(np), &
            p_orig(np))
    endif
    if (.not.allocated(p_saved)) then
        allocate(p_saved(np,noptfuncstore))
    endif

    !Set up positivep array (determines those parameters which must remain positive)
    call setuppositivep

    !Store initial values of parameters so they can be returned to their
    !original units once CG optimization is complete
    do i=1,np
        p_orig(i)=p(i)
    enddo

    !Determine the number of parameters to be optimized
    nfreep=0
    do i=1,np
        if ((freep(i).eq.1).or.(freep(i).eq.2)) then
            nfreep=nfreep+1
        endif
    enddo

    if (.not.allocated(popt)) then
        allocate(popt(nfreep))
    endif
    if (.not.allocated(d_smsno)) then
        allocate(d_smsno(nfreep))
        allocate(urp(nfreep,3))
    endif

    !Set up 'popt', the array containing the (rescaled) potential parameters for
    !use in the CG optimizer
    ip=0
    do i=1,np
        if ((freep(i).eq.1).or.(freep(i).eq.2)) then
            ip=ip+1
            if (p_orig(i).ne.0d0) then
                popt(ip)=p(i)*(poptStrt/p_orig(i))
            else
                popt(ip)=p(i)
            endif
        endif
    enddo

    !---- Set up variables necessary for call to CG optimizer ----
    lv=71+nfreep*(nfreep+15)/2 + 500
    if (.not.allocated(vv)) then
      allocate(vv(lv))
    endif
    d_smsno=1d0
    call deflt(2, iv, liv, lv, vv) !Initializes many variables automatically
    iv(1)=12 !fresh start; already read defaults (with changes below)
    iv(17)=maxfuncevals !max no. function evals
    iv(18)=10000 !max no. itns
    iv(20)=0 !don't print non default values
    iv(22)=0 !don't print final x and d
    iv(23)=0 !summary statistics
    iv(24)=0 !don't print initial x and d
    vv(32) = optdiff
    vv(42) = optfunc_err  !Relative error in calculating optfunc (important for gradient calc)
    opt_failed=.false.
    !-------------------------------------------------------------

    !---- CG optimizer ----
    if (nsteps.gt.1) then
       call smsno(nfreep,d_smsno,popt,Fnew,iv,liv,lv,vv,uip,urp,Fnew)
    endif
    !----------------------

    if (opt_failed.eqv..true.) then
        return
    endif
    lowestoptfuncrand=100000d0

    !Display optimization data
    if (nsteps.gt.1) then
       write(*,*) 'Optimization cycle completed.'
    endif
    call Fnew( nfreep, popt, tmpint, func, uip, urp, Fnew )
    print *,'Optimization function=',func
    print *,'Time:'
    call displaytime(6,hrs)

    !Restore parameters from popt() to the p() array
    ip=0
    do i=1,np
        if ((freep(i).eq.1).or.(freep(i).eq.2)) then
            ip=ip+1
            if (p_orig(i).ne.0d0) then
                p(i)=popt(ip)*(p_orig(i)/poptStrt)
            else
                p(i)=popt(ip)
            endif
        endif
    enddo
    !Ensure meam parameters have the correct signs and are in the correct range
    call parachange1(signflipped,cutoff,nf)

    !Old code: retain for re-inclusion at later date:
    !    call optimizestructure(.false.) !First index: optimize all structures?
    !    call initializestruc(0)
    !    call setupnntables
    !    !Relax all structures using the potential (and not just the ones
    !    !contributing to the optimization function)
    !    call optimizestructure(.true.)
    !    call initializestruc(0)
    !    call setupnntables

    !Save potential and record energies and forces
    if (func.ne.0d0) then !Check optimization function properly computed

        !Re-order 'noptfuncstore' best optfuncs
        if (n_optfunc.gt.1) then
            largest=.true. !Keep track of whether current optfunc is larger than
                           !the stored values
            do i=1,n_optfunc-1
                if ((func.lt.bestoptfuncs(i)).and.(largest.eqv..true.)) then
                    if (i.eq.1) then
                        print *,'New best optfunc !!!'
                    endif
                    !Move all the optfunc values greater than the new value
                    !up by one element, and do the same for the stored MEAM
                    !parameters
                    do j=min(n_optfunc-1,noptfuncstore-1),i,-1
                        bestoptfuncs(j+1)=bestoptfuncs(j)
                        timeoptfunc(j+1)=timeoptfunc(j)
                        do k=1,np
                            p_saved(k,j+1)=p_saved(k,j)
                        enddo
                    enddo
                    bestoptfuncs(i)=func
                    timeoptfunc(i)=hrs
                    !Record the potential parametets
                    do k=1,np
                        p_saved(k,i)=p(k)
                    enddo
                    largest=.false.
                endif
            enddo
            if (largest.eqv..true.) then
                !Add the potential to the back of the stored potentials
                if (n_optfunc.le.noptfuncstore) then
                    bestoptfuncs(n_optfunc)=func
                    timeoptfunc(n_optfunc)=hrs
                    !Record the potential parameters
                    do k=1,np
                        p_saved(k,n_optfunc)=p(k)
                    enddo
                endif
            endif
        else
            if (n_optfunc.le.noptfuncstore) then
                !Save optfunc as the best optfunc
                bestoptfuncs(1)=func
                timeoptfunc(1)=hrs
                !Save MEAM parameters as the best MEAM parameters
                do i=1,np
                    p_saved(i,1)=p(i)
                enddo
                largest=.false.
            endif
        endif
        if (n_optfunc.le.noptfuncstore) then
            n_optfunc=n_optfunc+1
        endif

        !Save best optfuncs to file 'bestoptfuncs' and write to standard output
        print *
        if (n_optfunc.gt.2) then
            open(50,file='bestoptfuncs')
            if (n_optfunc.gt.noptfuncstore) then
                print *,'Top ',noptfuncstore,' optfuncs:'
                write(50,*) 'Top ',noptfuncstore,' optfuncs:'
            else
                print *,'Best optfuncs so far:'
            endif
            do i=1,n_optfunc-1
                print *,'   ',i,':',bestoptfuncs(i),' time: ',INT(timeoptfunc(i)),' hours'
                write(50,*) '   ',i,':',bestoptfuncs(i),' time:',INT(timeoptfunc(i)),' hours'
            enddo
            print *
            write(50,*)
            call displaytime(50,hrs)
            close(50)
        endif

        !Save best MEAM paras to files
        do i=1,n_optfunc-1
            do j=1,np
                p(j)=p_saved(j,i)
            enddo
            !Generate a filename of the form 'potparas_best1', etc
            filename="potparas_best"
            if (i.lt.10) then
                write(string1,'(I1)') i
            elseif (i.lt.100) then
                write(string1,'(I2)') i
            else
                print *,'ERROR: more than 100 files set to save; code needs'
                print *,'changing, STOPPING.'
                stop
            endif
            filename=trim(filename)//trim(string1)
            !Write potential parameters to file
            open(2,file=trim(adjustl(filename)))
            call writemeamp
            close(2)
 
        enddo

        !Save fit data and true data for best MEAM paras to files and
        !plot potential functions
        do i=1,n_optfunc-1

            !Copy saved configurations for ith optfunc to popt
            ip=0
            do j=1,np
                if ((freep(j).eq.1).or.(freep(j).eq.2)) then
                    ip=ip+1
                    if (p_orig(j).ne.0d0) then
                        popt(ip)=p_saved(j,i)*(poptStrt/p_orig(j))
                    else
                        popt(ip)=p(j)
                    endif
                endif
            enddo

            !Generate a filename of the form 'potparas_best1', etc
            filename="datapnts_best"
            if (i.lt.10) then
                write(string1,'(I1)') i
            elseif (i.lt.100) then
                write(string1,'(I2)') i
            else
                print *,'ERROR: more than 100 files set to save; code needs'
                print *,'changing, STOPPING.'
                stop
            endif
            filename=trim(filename)//trim(string1)
            open(2,file=trim(adjustl(filename)))

            !Calculate and then write the fitted datapoints to the above file
            call Fnew( nfreep, popt, tmpint, func, uip, urp, Fnew )
            call writemeamdatapoints(2,.false.)
            close(2)

            !Plot potential functions
            if (lmax.eq.0) then
               !Setup name of LAMMPS output file
               filename=""
               do j=1,maxspecies
                   write(string1,'(A2)') element(speciestoZ(j))
                   filename=trim(filename)//trim(string1)
               enddo
               string2='.eam.alloy_'
               if (i.lt.10) then
                   write(string3,'(I1)') i
               elseif (i.lt.100) then
                   write(string3,'(I2)') i
               else
                   print *,'ERROR: more than 100 files set to save; code needs'
                   print *,'changing, STOPPING.'
                   stop
               endif
               filename=trim(filename)//trim(string2)//trim(string3)
               open(58,file=trim(adjustl(filename)))
            endif
            !Calculate the largest encountered background density, necessary to
            !set the upper range on the output embedding functions
            maxlargestrho=0d0
            do j=1,nstruct
                call getlargestbkgrnddens
                maxlargestrho=MAX(maxlargestrho,largestrho)
            enddo
            !Plot the potential and write LAMMPS output
            call plotfunctions(.false.)
            if (lmax.eq.0) then
               close(58)
            endif
        enddo

        !Check if optimization termination conditions have been met
        if (n_optfunc.gt.2) then
            if ( (n_optfunc.gt.noptfuncstore).and.(bestoptfuncs(noptfuncstore)- &
                  bestoptfuncs(1)).lt.optAcc  ) then
               print *,'Best and worst optimization functions within ',optAcc
               print *
               print *,' ::: Optimization completed :::'
               print *
               stop
           !else
           !   print *,'bestoptfuncs(noptfuncstore)- bestoptfuncs(1)=', &
           !            bestoptfuncs(noptfuncstore)- bestoptfuncs(1),' > ', &
           !           optAcc,', so continuing optimization'
            endif
         endif
         if ( hrs.ge.stopTime ) then
            print *,'Number of hours elapsed: ',hrs,' greater than STOPTIME=',stopTime
            print *
            print *,' ::: Optimization completed :::'
            print *
            stop
         endif

    endif

    !Reset popt() array to best performing potential
    ip=0
    do j=1,np
        if ((freep(j).eq.1).or.(freep(j).eq.2)) then
            ip=ip+1
            if (p_orig(j).ne.0d0) then
                popt(ip)=p_saved(j,1)*(poptStrt/p_orig(j))
            else
                popt(ip)=p(j)
            endif
        endif
    enddo

    !Check if 'smallest optfunc' criteria for convergence met
    if (func.gt.maxoptfuncallowed) then
        write(*,'(A9,f10.5,A2,f10.5,A15)') &
            ' Optfunc=',func,' >',maxoptfuncallowed, &
            '. Trying again.'
        lowestoptfuncrand=100000d0
        deallocate(popt)
        !call resetxyz
        return
    endif

    !...if so, stop and write out data to screen
    if (func.ne.0d0) then
        nsteps=1
        !With nsteps set to one, the F subroutine prints out the final
        !data and stops the program.
        call Fnew( nfreep, popt, tmpint, func, uip, urp, Fnew )
        deallocate(popt)
    endif

end subroutine optimizeparameters
subroutine optimizestructure(allstrucs)

    !--------------------------------------------------------------c
    !
    !     Relax the atomic coordinates of the first POSCAR file
    !     in 'fitdbse'. allstrucs determines if only the files
    !     which have been specified to relax and which have weights
    !     >0 are relaxed (.false.) or whether all files which have
    !     been specified to relax are relaxed regardless of the
    !     weight.
    !
    !     Called by:     program MEAMfit
    !     Calls:         initializestruc, initializemeam,
    !                 setuppositivep,parachange1,meamforce,
    !                 broyden
    !     Arguments:     allstrucs
    !     Returns:       -- nothing --
    !     Files read:    poscarfiles
    !     Files written: crystalstrucfile
    !
    !     Marcel Sluiter and Andy Duff, March 14 2011
    !
    !--------------------------------------------------------------c


    use m_geometry
    use m_datapoints
    use m_optimization
    use m_poscar
    use m_generalinfo
    use m_filenames
    use m_meamparameters

    implicit none
    logical allstrucs,firstrlx
    character(len=1) :: back
    integer, parameter:: maxiter=50!1000
    real(8), parameter:: forcetolerance=1d-3 !eV/angstrom
    real(8), parameter:: maxmovement=0.5, & !angstrom. Try bigger step
        !for Fe: 1.0, recently 1
    assumedspringconstant=1d0
    logical, allocatable:: signflipped(:)
    integer i,k,j,iter,nstor,nfreep,nf
    real(8) force(3)
    real(8) mixparam,latticeconst,poscarcell(9),magforce,maxmagforce, &
        aux,dx(3),dx2,scaledispl
    real(8), allocatable:: x(:),work(:),f(:),cutoff(:)
    character*80 poscarheader,poscarspecies

    print *,'Optimizing structure'

    mixparam=0.001d0
    do istr=1,nstruct
        optimizeforce_backup(istr)=optimizeforce(istr)
    enddo

    optimizeforce=0
    do istr=1,nstruct
        if (rlxstruc(istr).eq.1) then
            if ((allstrucs.eqv..true.).or.(weights(istr).gt.0d0)) then
                optimizeforce(istr)=1
            endif
        endif
    enddo

    call setupnntables
    firstrlx=.true.
    do istr=1,nstruct
        if (rlxstruc(istr).eq.1) then

            if ( (weights(istr).gt.0d0).or. (allstrucs.eqv..true.) ) then

                if (firstrlx) then
                    firstrlx=.false.
                    write(*,*) ' Relaxing all structures.'
                endif

                write (*,'(A10,A)') 'Relaxing ',strucnames(istr)

                !Update neighborhood
                call initializestruc(istr)

                if (gn_forces(istr).gt.1000) then
                    print *,'ERROR: number of atomic forces in first ', &
                        'POSCAR'
                    print *,'file is larger than 1000 (=',gn_forces(istr), &
                        ' Need to increase ', &
                        'size of force_store array in ', &
                        'optimizestructure subroutine. STOPPING.'
                    stop
                endif
                nstor = min(5,3*gn_forces(istr),maxiter)
                !               print *,'allocating x of size:',3*gn_forces(istr)
                if (allocated(x)) deallocate(x)
                if (allocated(f)) deallocate(f)
                if (allocated(work)) deallocate(work)

                allocate(x(3*gn_forces(istr)),f(3*gn_forces(istr)), &
                    work(3*gn_forces(istr)+ &
                    2*nstor*(1+nstor+3*gn_forces(istr))) )

                !Set up meam parameters
                if (allocated(p)) deallocate(p)
                if (allocated(scaledstep)) deallocate(scaledstep)
                if (allocated(failedtries)) deallocate(failedtries)
                if (allocated(positivep)) deallocate(positivep)
                if (allocated(cutoff)) deallocate(cutoff)
                if (allocated(signflipped)) deallocate(signflipped)
                allocate(p(np),scaledstep(np),failedtries(np), &
                    positivep(np),cutoff(np),signflipped(np))
                call setuppositivep
                call parachange1(signflipped,cutoff,nf) !Ensure meam
                !parameters have the correct signs and are in
                !the correct range
                ! call setupsplines !Prepare splines for pair-potentials

                iter=1  !Loop over structural optimization steps
                do

                    !Calculate forces
                    maxmagforce=0d0 !Largest magnitude of force on any
                    !atom
                    do j=1,gn_forces(istr)
                        !                  print *,'j=',j
                        call meamforce(istr,j,force) !Calculate the
                        !force acting on atom j in structure i.
                        magforce= &
                            sqrt( force(1)**2 + force(2)**2 + force(3)**2 )
                        !                        print *,'magforce=',magforce
                        maxmagforce=max(magforce,maxmagforce)
                        x(3*(j-1)+1)=gxyz(1,j,istr)
                        x(3*(j-1)+2)=gxyz(2,j,istr)
                        x(3*(j-1)+3)=gxyz(3,j,istr)
                        f(3*(j-1)+1)=force(1)
                        f(3*(j-1)+2)=force(2)
                        f(3*(j-1)+3)=force(3)
                    enddo
                    !                  print *,'done force calc'
                    back=char(8)
                    if (iter.eq.1) then
                        write(*,'(A30,F15.10,A10,F15.10,A7)') &
                            '  Maximum magnitude of force:', &
                            maxmagforce,' eV/Ang > ',forcetolerance, &
                            ' eV/Ang'
                    else
                        write(*,'(256A1)') (back, k=1, 89)
                        write(*,'(A30,F15.10,A10,F15.10,A7)') &
                            '  Maximum magnitude of force:', &
                            maxmagforce,' eV/Ang > ',forcetolerance, &
                            ' eV/Ang'
                    endif

                    !Determine good mixparam for 1st structural
                    !optimization only
                    if (iter.eq.1) then
                        !                    aux=min(maxmovement,maxmagforce/
                        !     +assumedspringconstant)
                        !                    mixparam=aux/maxmagforce
                        !                    print *,mixparam
                        !                    stop
                        mixparam=0.05d0!0.04 gives the lowest second value
                    endif

                    call broyden(iter,3*gn_forces(istr), &
                        nstor,mixparam,work,x,f)

                    !Check that position update does not move atoms more
                    !than maxmovement
                    scaledispl = 1d0
                    aux = 0d0
                    do i=1,gn_forces(istr)
                        dx(1) = gxyz(1,i,istr)-x(3*(i-1)+1)!f(3*(i-1)+1)
                        dx(2) = gxyz(2,i,istr)-x(3*(i-1)+2)!f(3*(i-1)+2)
                        dx(3) = gxyz(3,i,istr)-x(3*(i-1)+3)!f(3*(i-1)+3)
                        dx2   = sqrt( dot_product(dx,dx) )
                        aux   = max( aux, dx2 )
                        if (dx2.gt.maxmovement) then
                            !                         write(*,'(a,a,i3,a,f9.5)')
                            !     +                  '*WARNING* optimizestructure: movement too',
                            !     +                  ' great, atom:',i,' movement [Ang]=',dx2
                            scaledispl=min(scaledispl,maxmovement/dx2)
                        endif
                    enddo
                    if(abs(scaledispl-1d0).gt.1d-6) then
                        scaledispl=scaledispl*0.9d0 !arbitrary, but safe
                        !                    write(*,'(a,a,f9.5)')
                        !     +               'movement of atoms has been scaled ',
                        !     +               'by factor:',scaledispl
                    endif
                    do i=1,gn_forces(istr)
                        dx(1) =  gxyz(1,i,istr)-x(3*(i-1)+1)
                        dx(2) =  gxyz(2,i,istr)-x(3*(i-1)+2)
                        dx(3) =  gxyz(3,i,istr)-x(3*(i-1)+3)
                        dx=dx*scaledispl
                        gxyz(1,i,istr) = gxyz(1,i,istr) - dx(1) !update
                        gxyz(2,i,istr) = gxyz(2,i,istr) - dx(2) !positions
                        gxyz(3,i,istr) = gxyz(3,i,istr) - dx(3)
                    enddo

                    !Update the POSCAR file with the new atomic positions.
                    !This is necessary at the moment, because the
                    !subroutine 'initializestruc' reads in from the POSCAR
                    !files, and this subroutine has to be called to reset
                    !the atomic environments (neighbors) each time the
                    !atoms are moved. Note: update so that this is not
                    !necessary.
                    open(unit=1,file=trim(poscarfiles(istr)))
                    read(1,'(a)') poscarheader
                    read(1,*) latticeconst
                    read(1,*) poscarcell
                    read(1,'(a)') poscarspecies
                    rewind(unit=1)
                    write(1,'(a)') trim(poscarheader)
                    write(1,*) latticeconst
                    write(1,'(3f12.6)') poscarcell
                    write(1,'(a)') trim(poscarspecies)
                    write(1,'(a)') 'Cartesian'
                    do j=1,gn_forces(istr)
                        write(1,*) gxyz(1:3,j,istr)/latticeconst
                        !                     print *,gxyz(1:3,j,1)/latticeconst
                    enddo
                    close(1)
                    !                  stop

                    iter=iter+1
                    if ((iter.gt.maxiter).or. &
                        (maxmagforce.lt.forcetolerance)) then
                        exit
                    endif

                    !Update neighborhood
                    call initializestruc(istr)

                enddo !end loop over structural optimization
                if (maxmagforce.ge.forcetolerance) then
                    write(*,'(256A1)') (back, k=1, 89)
                    write(*,'(A30,F15.10,A10,F15.10,A7)') &
                        '  Maximum magnitude of force:', &
                        maxmagforce,' eV/Ang > ',forcetolerance, &
                        ' eV/Ang'
                    print *
                    write(*,'(A20,A14,I4,A1,I4,A1)') '  Structural optim.', &
                        'failed (iter:',maxiter,'/',maxiter,')'
                else
                    write(*,'(256A1)') (back, k=1, 89)
                    write(*,'(A30,F15.10,A10,F15.10,A7)') &
                        '  Maximum magnitude of force:', &
                        maxmagforce,' eV/Ang < ',forcetolerance, &
                        ' eV/Ang'
                    print *
                    write(*,'(A20,A12,I4,A1,I4,A1)') '  Structural optim.', &
                        'done (iter:',iter,'/',maxiter,')'
                endif

                !               deallocate(x,f,work,p,scaledstep,failedtries,
                !     +                    positivep,cutoff,signflipped)

            endif

        endif

    enddo

    !De-allocate relevant arrays (ANDY should you not put
    !this
    !stuff outside the loop???)
    !      if (allocated(zz)) deallocate(zz)
    !      if (allocated(gn_inequivalentsites))
    !     +     deallocate(gn_inequivalentsites)
    !      if (allocated(gspecies)) deallocate(gspecies)
    !      IF (allocated(gn_neighbors)) deallocate(gn_neighbors)
    !      if (allocated(gneighborlist)) deallocate(gneighborlist)
    !      if (allocated(gxyz)) deallocate(gxyz)
    !      if (allocated(gn_forces)) deallocate(gn_forces)
    !      if (allocated(gn_C)) deallocate(gn_C)
    !      if (allocated(optforce)) deallocate(optforce)

    do istr=1,nstruct
        optimizeforce(istr)=optimizeforce_backup(istr)
    enddo

end subroutine optimizestructure
subroutine optionalarg(line,nargs)

    !--------------------------------------------------------------c
    !
    !     Parses the string 'line', searches for number of strings
    !     separated by spaces or commas.
    !
    !     Andrew Ian Duff (Mar 25 2013): words within '' are only
    !     counted as one. This was important, because when the words
    !     of line are read into individual words with the
    !     read(line,*) words(1:nargs) command in the
    !     readoutcarenergy subroutine, words within
    !     apostrophes are counted as a single word, irrespective
    !     of any spaces contained within. This resulted in the above
    !     command trying to read in a non-existent word from
    !     line into the word array and thus crashing
    !
    !     Called by:     readposcar,readoutcarenergy,
    !                 readoutcarforces
    !     Calls:         -
    !     Arguments:     line
    !     Returns:       nargs
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter, Dec 28 2005,
    !     Edited by Andrew Ian Duff, Mar 25 2013.
    !
    !--------------------------------------------------------------c

    implicit none

    integer nargs,length,i
    character line*(*)
    logical newarg,donotcount

    intent(in) line
    intent(out) nargs

    !Analyze LINE
    length=len_trim(line)
    !Count number of sections in LINE that are separated by commas or
    !spaces
    nargs=0
    newarg=.true.
    donotcount=.false. !this flag gets switched on when ' is found.
    !No 'sections' are counted until another ' is
    !found.
    do i=1,length
        if (line(i:i).eq.CHAR(39)) then
            if (donotcount.eqv..false.) then
                donotcount=.true.
            else
                donotcount=.false.
            endif
        endif
        if (donotcount.eqv..false.) then
            if ((line(i:i).eq.',').or.(line(i:i).eq.' ').or. &
                (line(i:i).eq.CHAR(9))) then
                newarg=.true.
            else
                if(newarg) nargs=nargs+1
                newarg=.false.
            endif
        endif
    enddo

end subroutine optionalarg
function pairpotAcklandFeFe(distance)

    implicit none

    real(8) distance,pairpotAcklandFeFe,cutoff

    if (distance.lt.1d0) then
        pairpotAcklandFeFe=( 9734.2365892908d0 / distance ) * &
            ( 0.18180d0*exp(-2.8616724320005d1*distance) & !E+01
            + 0.50990d0*exp(-8.4267310396064d0*distance) &
            + 0.28020d0*exp(-3.6030244464156d0*distance) &
            + 0.02817d0*exp(-1.8028536321603d0*distance) )
    else
        pairpotAcklandFeFe=exp( 7.4122709384068d0 - &
            0.64180690713367d0*distance - &
            2.6043547961722d0*(distance**2) + &
            0.6262539393123d0*(distance**3) ) * &
            cutoff(distance,1d0,2.05d0) &
            - 27.444805994228d0 * ((2.2d0 - distance) **3) * &
            cutoff(distance,2.05d0,2.2d0) &
            + 15.738054058489d0 * ((2.3d0 - distance) **3) * &
            cutoff(distance,2.05d0,2.3d0) &
            + 2.2077118733936d0 * ((2.4d0 - distance) **3) * &
            cutoff(distance,2.05d0,2.4d0) &
            - 2.4989799053251d0 * ((2.5d0 - distance) **3) * &
            cutoff(distance,2.05d0,2.5d0) &
            + 4.2099676494795d0 * ((2.6d0 - distance) **3) * &
            cutoff(distance,2.05d0,2.6d0) &
            - 0.77361294129713d0 * ((2.7d0 - distance) **3) * &
            cutoff(distance,2.05d0,2.7d0) &
            + 0.80656414937789d0  * ((2.8d0 - distance) **3) * &
            cutoff(distance,2.05d0,2.8d0) &
            - 2.3194358924605d0  * ((3.d0 - distance) **3) * &
            cutoff(distance,2.05d0,3.d0) &
            + 2.6577406128280d0  * ((3.3d0 - distance) **3) * &
            cutoff(distance,2.05d0,3.3d0) &
            - 1.0260416933564d0  * ((3.7d0 - distance) **3) * &
            cutoff(distance,2.05d0,3.7d0) &
            + 0.35018615891957d0  * ((4.2d0 - distance) **3) * &
            cutoff(distance,2.05d0,4.2d0) &
            - 0.058531821042271d0  * ((4.7d0 - distance) **3)* &
            cutoff(distance,2.05d0,4.7d0) &
            - 0.0030458824556234d0  * &
            ((5.3d0 - distance) **3)* &
            cutoff(distance,2.05d0,5.3d0)
    endif

end function pairpotAcklandFeFe
subroutine pairpotential(species1,species2,distance,pairpot)

    !------------------------------------------------------------------c
    !
    !     Calculates the pairpotential, 'pairpot', for two atoms
    !     of species, 'species1' and 'species2' which are a
    !     distance, 'distance', apart from one another.
    !
    !     Called by:     meamenergy
    !     Calls:         -
    !     Arguments:     species1,species2,distance,pairpotparameter
    !     Returns:       pairpot
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter, Feb 9 2006
    !
    !------------------------------------------------------------------c

    use m_atomproperties
    use m_meamparameters
    use m_filenames
    use m_generalinfo
    use m_optimization

    implicit none

    integer species1,species2,i,index,isp,jsp!,splnPosn,
    real(8) distance,pairpot,tmp,rsBier,x,eps,B0,B1,B2,B3,k1,k2,k3,k4,depsdr, &
            rdiffcb,rBiersq,rBiercb,cutoffMinsq,cutoffMincb!,diff
    real(8), parameter :: rBier=0.9d0
    real(8), parameter :: twothrds=0.666666666666667d0
    real(8), parameter :: e2ov4pieps0=14.3977847299112d0 !e^2 / 4 pi eps_0, in eV/Angstrom
    ! ( = (1.602*10^-19)^2 / (4*3.142*8.854*10^-12) * (6.241 * 10^18) * 10^10
    !                                             1 J = 6.241*10^18 eV   1 m = 10^10 Angstron
    intent(in) species1,species2!distance
    intent(out) pairpot
    !print *,'species1=',species1,', species2=',species2
    !print *,'Z1=',speciestoZ(species1),' Z2=',speciestoZ(species2)

    if (distance.lt.rBier) then

      rsBier = 0.88534d0 * 0.52917721092d0 / sqrt( ( speciestoZ(species1) )**twothrds &
                         + ( speciestoZ(species2) )**twothrds  )
      x = distance/rsBier
      eps = 0.1818d0*exp(-3.2d0*x)  + 0.5099d0*exp(-0.9423d0*x) + &
            0.2802d0*exp(-0.4029*x) + 0.02817d0*exp(-0.2016*x)
      pairpot = e2ov4pieps0 * speciestoZ(species1) * speciestoZ(species2) * eps / distance
    elseif (distance.lt.cutoffMin) then

      !Here we interpolate between the lower radius Biersack result and the
      !larger radius sum over cubic terms. We must first evaluate the pairpotential 
      !and its derivative at the lower and upper radial values.

      !Evaluate k1: value of pairpotential at R=rBier
      rsBier = 0.88534d0 * 0.52917721092d0 / sqrt( ( speciestoZ(species1))**twothrds &
                         + ( speciestoZ(species2) )**twothrds  )
      x = rBier/rsBier
      eps = 0.1818d0*exp(-3.2d0*x)  + 0.5099d0*exp(-0.9423d0*x) + &
            0.2802d0*exp(-0.4029d0*x) + 0.02817d0*exp(-0.2016d0*x)
      k1 = e2ov4pieps0 * speciestoZ(species1) * speciestoZ(species2) * eps / rBier

      !Evaluate k2: value of d(pairpotential)/dr at R=rBier
      depsdr = (1d0/rsBier)*( 0.1818d0*(-3.2d0)*exp(-3.2d0*x) + 0.5099d0*(-0.9423d0)*exp(-0.9423d0*x) + &
            0.2802d0*(-0.4029d0)*exp(-0.4029d0*x) + 0.02817d0*(-0.2016d0)*exp(-0.2016d0*x) )

      k2 = - e2ov4pieps0 * speciestoZ(species1) * speciestoZ(species2) * eps / rBier**2 + &
           e2ov4pieps0 * speciestoZ(species1) * speciestoZ(species2) * depsdr / rBier

      !Evaluate k3: value of pairpotential at R=cutoffMin
      k3 = 0
      if (lookuptables.eqv..true.) then
         print *,'Need to update code in pairpotential subroutine, stopping.'
         stop
      endif
     
      do isp=1,maxspecies
         do jsp=isp,maxspecies
            if ( ((species1.eq.isp).and.(species2.eq.jsp)).or. &
                 ((species1.eq.jsp).and.(species2.eq.isp)) ) then
               if (typepairpot(isp,jsp).eq.2) then
                  do i=1,31,2
                     if (cutoffMin.le.pairpotparameter(i+1,isp,jsp)) then
                         k3=k3 + pairpotparameter(i,isp,jsp)  * &
                             ((pairpotparameter(i+1,isp,jsp)   - cutoffMin) **3)
                     endif
                  enddo
               else
                  print *,'typepairpot(',isp,',',jsp,')=',typepairpot(isp,jsp), &
                          'not implemented, stopping.'
                  stop
               endif
            endif
         enddo
      enddo

      !Evaluate k4: value of d(pairpotential)/dr at R=cutoffMin
      k4 = 0
      do isp=1,maxspecies
         do jsp=isp,maxspecies
            if ( ((species1.eq.isp).and.(species2.eq.jsp)).or. &
                 ((species1.eq.jsp).and.(species2.eq.isp)) ) then
               do i=1,31,2
                  if (cutoffMin.le.pairpotparameter(i+1,isp,jsp)) then
                     k4=k4 - 3d0*pairpotparameter(i,isp,jsp) * &
                         ((pairpotparameter(i+1,isp,jsp)-cutoffMin)**2)
                  endif
               enddo
            endif
         enddo
      enddo

      !Now calculate B0-B3 from the values k1-k4
      !...but first, preliminary quantities:      
      rdiffcb=(rBier-cutoffMin)**3
      rBiersq=rBier**2
      rBiercb=rBier**3
      cutoffMinsq=cutoffMin**2
      cutoffMincb=cutoffMin**3
     !print *,'rBier=',rBier
     !print *,'rdiffcb=',rdiffcb,', rBiersq=',rBiersq
     !print *,'rBiercb=',rBiercb,', cutoffMin=',cutoffMin
     !print *,'cutoffMinqs=',cutoffMinsq,', cutoffMincb=',cutoffMincb
      B0=(-1d0/rdiffcb) * ( -k3*rBiercb + &
            3d0*k3*rBiersq*cutoffMin + &
            k4*rBiercb*cutoffMin - &
            3d0*k1*rBier*cutoffMinsq + &
            k2*rBiersq*cutoffMinsq - &
            k4*rBiersq*cutoffMinsq + &
            k1*cutoffMincb - k2*rBier*cutoffMincb )
      B1=(-1d0/rdiffcb) * ( -k4*rBiercb + &
            6d0*k1*rBier*cutoffMin - &
            6d0*k3*rBier*cutoffMin - &
            2d0*k2*rBiersq*cutoffMin - &
            k4*rBiersq*cutoffMin + &
            k2*rBier*cutoffMinsq + 2d0*k4*rBier*cutoffMinsq + &
            k2*cutoffMincb )
      B2=(-1d0/rdiffcb) * ( -3d0*k1*rBier + 3d0*k3*rBier + &
            k2*rBiersq + 2d0*k4*rBiersq - &
            3d0*k1*cutoffMin + 3d0*k3*cutoffMin + &
            k2*rBier*cutoffMin - k4*rBier*cutoffMin - &
            2d0*k2*cutoffMinsq - k4*cutoffMinsq )
      B3=-(2d0*k1-2d0*k3-k2*rBier-k4*rBier+k2*cutoffMin+k4*cutoffMin) / &
          (rdiffcb)
     !print *,'B0=',B0
     !print *,'B1=',B1
     !print *,'B2=',B2
     !print *,'B3=',B3
      pairpot = B0 + B1*distance + B2*(distance**2) + B3*(distance**3)

     !distance=rBier
     !pairpot = B0 + B1*distance + B2*(distance**2) + B3*(distance**3)
     !print *,'interpolated value at r=rBier, =',pairpot
     !print *,'... should equal k1=',k1
     !print *,'interpolated derivative at r=rBier=',B1+2d0*B2*distance+3d0*B3*(distance**2)
     !tmp=pairpot
     !distance=rBier+0.0001d0
     !pairpot = B0 + B1*distance + B2*(distance**2) + B3*(distance**3)
     !print *,'by finite difference=',(pairpot-tmp)/0.0001d0
     !print *,'... should equail k2=',k2

     !distance=cutoffMin
     !pairpot = B0 + B1*distance + B2*(distance**2) + B3*(distance**3)
     !print *,'interpolated value at r=cutoffMin, =',pairpot
     !print *,'... should equal k3=',k3
     !print *,'interpolated derivative at r=rBier=',B1+2d0*B2*distance+3d0*B3*(distance**2)
     !tmp=pairpot
     !distance=cutoffMin+0.0001d0
     !pairpot = B0 + B1*distance + B2*(distance**2) + B3*(distance**3)
     !print *,'by finite difference=',(pairpot-tmp)/0.0001d0
     !print *,'... should equail k4=',k4
     !stop

    else

!      if (lookuptables.eqv..true.) then
!         !Add code for reading in data direct from tables here.
!      else

      do isp=1,maxspecies
         do jsp=isp,maxspecies
            if ( ((species1.eq.isp).and.(species2.eq.jsp)).or. &
                 ((species1.eq.jsp).and.(species2.eq.isp)) ) then
               if (typepairpot(isp,jsp).eq.2) then
                  pairpot=0d0
                  do i=1,31,2
                     if (distance.le.pairpotparameter(i+1,isp,jsp)) then
                         pairpot=pairpot + pairpotparameter(i,isp,jsp)  * &
                             ((pairpotparameter(i+1,isp,jsp)   - distance) **3)
                     endif
                  enddo
               endif
            endif
         enddo
      enddo

!      endif

    endif

end subroutine pairpotential



!     Old code relating to Ackland Fe-Fe potential:
!
!
!
!                 !Ackland Fe-Fe pair-potential
!
!                 !              pairpot=pairpotAcklandFeFe(distance)
!                 !              tmp=pairpot
!                 tmp=((distance/p_rmax)**0.5)*dble(narrpairpotXX)
!                 index=tmp
!                 call splint(r_pairpotXXarr,pairpotXXarr, &
!                     secderpairpotXX,narrpairpotXX,narrpairpotXX, &
!                     index,distance,pairpot)
!                 !              print *,pairpot
!                 !              r_pairpotXXarr(i) =
!                 !     +         ( (dble(i)/dble(narrpairpotXX)) )**2 *
!                 !     +          p_rmax
!                 !              tmp=((distance/p_rmax)**0.5)*dble(narrpairpotXX)
!                 !              i=tmp
!                 !              print *,pairpotXXarr(i),pairpotXXarr(i+1)
!                 !              stop
!
!                 !              write(91,*) distance,tmp,pairpot
!
!                 !            if (distance.lt.1d0) then
!                 !              pairpot=( 9734.2365892908d0 / distance ) *
!                 !    +               ( 0.18180d0*exp(-2.8616724320005d1*distance) !E+01
!                 !    +               + 0.50990d0*exp(-8.4267310396064d0*distance)
!                 !    +               + 0.28020d0*exp(-3.6030244464156d0*distance)
!                 !    +               + 0.02817d0*exp(-1.8028536321603d0*distance) )
!                 !            else
!                 !              pairpot=exp( 7.4122709384068d0 -
!                 !    +               0.64180690713367d0*distance -
!                 !    +               2.6043547961722d0*(distance**2) +
!                 !    +               0.6262539393123d0*(distance**3) ) *
!                 !    +               cutoff(distance,1d0,2.05d0)
!                 !    +              - 27.444805994228d0 * ((2.2d0 - distance) **3) *
!                 !    +               cutoff(distance,2.05d0,2.2d0)
!                 !    +              + 15.738054058489d0 * ((2.3d0 - distance) **3) *
!                 !    +               cutoff(distance,2.05d0,2.3d0)
!                 !    +              + 2.2077118733936d0 * ((2.4d0 - distance) **3) *
!                 !    +               cutoff(distance,2.05d0,2.4d0)
!                 !    +              - 2.4989799053251d0 * ((2.5d0 - distance) **3) *
!                 !    +               cutoff(distance,2.05d0,2.5d0)
!                 !    +              + 4.2099676494795d0 * ((2.6d0 - distance) **3) *
!                 !    +               cutoff(distance,2.05d0,2.6d0)
!                 !    +              - 0.77361294129713d0 * ((2.7d0 - distance) **3) *
!                 !    +                cutoff(distance,2.05d0,2.7d0)
!                 !    +              + 0.80656414937789d0  * ((2.8d0 - distance) **3) *
!                 !    +                 cutoff(distance,2.05d0,2.8d0)
!                 !    +              - 2.3194358924605d0  * ((3.d0 - distance) **3) *
!                 !    +                 cutoff(distance,2.05d0,3.d0)
!                 !    +              + 2.6577406128280d0  * ((3.3d0 - distance) **3) *
!                 !    +                 cutoff(distance,2.05d0,3.3d0)
!                 !    +              - 1.0260416933564d0  * ((3.7d0 - distance) **3) *
!                 !    +                 cutoff(distance,2.05d0,3.7d0)
!                 !    +              + 0.35018615891957d0  * ((4.2d0 - distance) **3) *
!                 !    +                 cutoff(distance,2.05d0,4.2d0)
!                 !    +              - 0.058531821042271d0  * ((4.7d0 - distance) **3)*
!                 !    +                 cutoff(distance,2.05d0,4.7d0)
!                 !    +              - 0.0030458824556234d0  *
!                 !    +              ((5.3d0 - distance) **3)*
!                 !    +                 cutoff(distance,2.05d0,5.3d0)
!                 !            endif
!                 !               print *,'pairpot analytic=',pairpot
!                 !               stop





!   Old code used when I tried using a spline rather than explicit calculation
!   of pairpotentials each time (did not give a speed increase so dropped this)
!
!
!
!
!           ! tmp=(splnNvals-1)*(distance-cutoffMin)/(cutoffMax-cutoffMin)
!           ! splnPosn=floor(tmp)+1 !Finds the lower bound array element
!           ! diff=tmp-floor(tmp)
!           ! pairpot=(1d0-diff)*pairpotStr(splnPosn,1,2) + &
!           !                diff*pairpotStr(splnPosn+1,1,2)
!
!              !if (pairpot.ne.0d0) then
!              !   write(81,*) 'distance=',distance
!              !   write(81,*) '   from array, pairpot(1,2)=',pairpot
!              !endif




!    Code for Fe-C hepburn potential:
!
!
!             !Following is Hepburn form, but HARD-WIRED
!             !               if (species1.eq.2) then
!             !                  print *,'distance=',distance
!             !               endif
!
!             !              do i=1,1000
!             !
!             !              distance=0.92d0+(dble(i-1)/999d0)*(4d0-0.92d0)
!
!             if (distance.lt.0.92d0) then
!                 print *,'need to implement full Fe-C pot'
!                 stop
!             elseif (distance.lt.1.7286d0) then
!                 pairpot=   1421.5171091881768d0 &
!                     -      3690.825251313584d0  *  distance &
!                     +      3710.621639067296    * (distance**2) &
!                     -      1679.6860367988738   * (distance**3) &
!                     +       285.39847077772725  * (distance**4)
!             elseif (distance.lt.1.88d0) then
!                 pairpot=   1521.0470658708155 &
!                     -      2389.4324313634993   *  distance &
!                     +      1252.1878589565072   * (distance**2) &
!                     -       218.9361400835967   * (distance**3)
!             elseif (distance.lt.2.25d0) then
!                 pairpot=    241.09303032117546 &
!                     -       346.952587401322    *  distance &
!                     +       165.7624100404569   * (distance**2) &
!                     -        26.307514389261673 * (distance**3)
!             elseif (distance.lt.2.42d0) then
!                 pairpot= -  581.6276060699361 &
!                     +       750.0082611201603   *  distance &
!                     -       321.77574485797965  * (distance**2) &
!                     +        45.9203604105067   * (distance**3)
!             elseif (distance.lt.2.7244d0) then
!                 pairpot=    364.0122288533755 &
!                     -       422.2725259748537   *  distance &
!                     +       162.63780352838967  * (distance**2) &
!                     -        20.803268843814134 * (distance**3)
!             elseif (distance.lt.3.1581d0) then
!                 pairpot= -  271.6958017654534 &
!                     +       277.7436580864025   *  distance &
!                     -        94.30544418165887  * (distance**2) &
!                     +        10.634019820362498 * (distance**3)
!             elseif (distance.lt.3.5d0) then
!                 pairpot=   4005.3336322295286 &
!                     -      4759.736033262177    *  distance &
!                     +      2117.9776780306693   * (distance**2) &
!                     -       418.29875898347865  * (distance**3) &
!                     +        30.94094273871914  * (distance**4)
!             else
!                 pairpot=0d0
!             endif
!
!
!
!    Hepburn for C-C:
!
!             !Following is Hepburn form, but HARD-WIRED
!             if (distance.lt.1d0) then
!                 print *,'tend to the pairpot subroutine...'
!                 stop
!             elseif (distance.lt.1.2857838598417968d0) then
!                 pairpot=   1199.5739471496045d0 &
!                     -      3835.3033425305593d0  *  distance &
!                     +      4786.499640264303d0  * (distance**2) &
!                     -      2705.687420612359d0 * (distance**3) &
!                     +       577.189423411637d0 * (distance**4)
!             elseif (distance.lt.1.8008513964923578d0) then
!                 pairpot=    286.8908260379106d0 &
!                     -       478.88759650017056d0 *  distance &
!                     +       267.62987770250356d0 * (distance**2) &
!                     -        49.910297272136546d0* (distance**3)
!             elseif (distance.lt.2.2863452818753887d0) then
!                 pairpot= -   16.399912881986573d0 &
!                     +        26.357988189333234d0 *  distance &
!                     -        12.929409795401792d0 * (distance**2) &
!                     +         2.020563142807547d0* (distance**3)
!             elseif (distance.lt.3.5d0) then
!                 pairpot=     11.60221676482102d0 &
!                     -        10.554683135163948d0 *  distance &
!                     +         3.3641528432894177d0 * (distance**2) &
!                     -         0.4199752489948345d0* (distance**3) &
!                     +         0.014225677158589412d0* (distance**4)
!             else
!                 pairpot=0d0
!             endif






! Old spline code:
!
!
!
!
!           ! tmp=(splnNvals-1)*(distance-cutoffMin)/(cutoffMax-cutoffMin)
!           ! splnPosn=floor(tmp)+1 !Finds the lower bound array element
!           ! diff=tmp-floor(tmp)
!           ! pairpot=(1d0-diff)*pairpotStr(splnPosn,1,2) + &
!           !                diff*pairpotStr(splnPosn+1,1,2)
!
!              !if (pairpot.ne.0d0) then
!              !   write(81,*) 'distance=',distance
!              !   write(81,*) '   from array, pairpot(1,2)=',pairpot
!              !endif

subroutine parachange1(signflipped,cutoff,nf)

    use m_optimization
    use m_generalinfo
    use m_meamparameters
    use m_atomproperties

    implicit none

    logical signflipped(np)
    integer i,ip,ispc,jspc,spc_cnt,ll,nf
    real(8) cutoff(np)

    !CutoffPenalty=0d0
    !CutoffPenCoeff=10000d0 !10 too small (order 10^-7 in optfunc in trial)

    !Check to see if any parameters ought to be positive, then
    !flip these parameters to be positive, but recording that the
    !sign has been flipped so that we can return the original negative
    !value to the conjugate-gradient minimizer.
    signflipped(1:np)=.false.
    do ip=1,np
        !         if (ip.eq.36) then
        !             print *,'for 36, positivep=',positivep(ip)
        !             print *,' and p(36)=',p(36)
        !c             stop
        !         endif
        if ((positivep(ip)).and.(p(ip).lt.0d0)) then
            signflipped(ip)=.true.
            p(ip)=abs(p(ip))
        endif
    enddo
    !Check to see if any of the radial cutoff parameters has exceeded
    !cutoffMax, or fallen below cutoffMin. If so, set these equal to 
    !cutoffMax/cutoffMin, but store the original value in
    !the cutoff array so the original value can be restored upon
    !returning the meam parameters to the conjugate gradient minimizer.
    cutoff(1:np)=0d0

    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    do i=2,12,2
                       ip=2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i
                       if (p(ip).gt.cutoffMax) then
                           !CutoffPenalty=CutoffPenalty+CutoffPenCoeff*((p(ip)-cutoffMax)**2)
                         !  nf=0 !SEE SMSNO in TOMS file (should be set to 0 if
                                !variables go out of bounds
                           cutoff(ip)=p(ip)
                           p(ip)=cutoffMax
                       endif
                       if ((p(ip).lt.cutoffMin).and.(p(ip).ne.0d0)) then
                           !CutoffPenalty=CutoffPenalty+CutoffPenCoeff*((cutoffMin-p(ip))**2)
                         !  nf=0
                           cutoff(ip)=p(ip)
                           p(ip)=cutoffMin
                       endif
                    enddo
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo
    
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                do i=2,32,2
                   ip=2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt
                   if (p(ip).gt.cutoffMax) then
                       !CutoffPenalty=CutoffPenalty+CutoffPenCoeff*((p(ip)-cutoffMax)**2)
                     !  nf=0
                       cutoff(ip)=p(ip)
                       p(ip)=cutoffMax
                   endif
                   if ((p(ip).lt.cutoffMin).and.(p(ip).ne.0d0)) then
                       !CutoffPenalty=CutoffPenalty+CutoffPenCoeff*((cutoffMin-p(ip))**2)
                     !  nf=0
                       cutoff(ip)=p(ip)
                       p(ip)=cutoffMin
                   endif
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

end subroutine parachange1
subroutine plotCutoffs

    !---------------------------------------------------------------c
    !
    !     Plots the cutoffs to the atomic density and pair- 
    !     potentials to files, so that they cen be plotted alongside
    !     the potential functions.
    !
    !     Called by:     -
    !     Calls:         -
    !     Arguments:     -
    !     Returns:       -
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff, 2014
    !
    !---------------------------------------------------------------c

    use m_filenames
    use m_meamparameters
    use m_optimization
    use m_atomproperties
    use m_generalinfo

    implicit none

    integer i,j,k,l,m

    j=1
    k=0
    l=1
    m=1
    do i=2*m3+lm1*m1+1,2*m3+lm1*m1+12*lm1*m2
        meamrhodecay(j,k,l,m)=p(i)
        if (j.lt.12) then
            j=j+1
        elseif (k.lt.lmax) then
            j=1
            k=k+1
        elseif (m.lt.maxspecies) then
            j=1
            k=0
            m=m+1
        elseif (l.lt.maxspecies) then
            j=1
            k=0
            l=l+1
            m=1
        endif
    enddo

    j=1
    if (maxspecies.eq.1) then
        do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,2*m3+(4+lm1)*m1+12*lm1*m2+32
            pairpotparameter(j,1,1)=p(i)
            j=j+1
        enddo
    else
        k=1
        l=1
        do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,m4
            pairpotparameter(j,k,l)=p(i)
            if (j.lt.32) then
                j=j+1
            elseif (l.lt.maxspecies) then
                j=1
                l=l+1
            elseif (k.lt.maxspecies) then
                j=1
                l=1
                k=k+1
            endif
        enddo
    endif

end subroutine plotCutoffs
subroutine plotfunctions(all)

    !     Currently for just one species. I have tested for a single
    !     species, namely for Fe, and got -4.013 eV for the bulk cell using
    !     the Fe.eam.alloy file as input to LAMMPS. Now I am generalising
    !     to FeC, and testing agreement between my code and LAMMPS.
    !     Currently doing for non FS-like radial densities (I.e., assume
    !     thi_YX=thi_XX and thi_YY=thi_XY). The energies will then be
    !     'wrong' using this code, but can at least test agreement with the
    !     LAMMPS values.
    !     NOTE: For the pair-potential, since the LAMMPS format demands the
    !     quantity V(r)*r @ r=0, and since V(r) must be evaluated first
    !     by 'pairpotential' before it is multiplied by r, we would get an
    !     infinity for this value. We could determine this value
    !     analytically, but then the code would not be general. Instead I
    !     interpolate the r=0 value from the first two lowest r values.
    !     Presumambly this would only introduce an error when atoms become
    !     infinatesimally close; I.e., never.
    !
    !     If all=.false., only produce the 58, file.

    use m_geometry
    use m_electrondensity
    use m_meamparameters
    use m_generalinfo
    use m_poscar
    use m_atomproperties
    use m_optimization
    use m_plotfiles
    use m_filenames

    implicit none

    logical all
    character*80 filename,string1,string2,string3
    integer i,j,jj,isp,jsp,Nr,Nrho,Npot,Nebf,gspeciesprev,l
    integer k,m,cnt,stringLen
    real(8) distmin,distmax,pairpot,sepn,drho,dr, &
        hartreeovereV,bohroverangstrom,deltar2,deltarho, &
        delta, LrgVal, &
        rho_out(5),meam_out(5),pairpot_out(5),sepn_out(5), &
        density_out(5)
    real(8) density(0:lmax)

    hartreeovereV=27.2116d0
    bohroverangstrom=0.529177249d0

    if (autorhomax.eqv..true.) then
       !Automatically choose rhomax based on the largest value of rho
       !encountered in the fit
        do i=-5,15
           if (maxlargestrho.lt.10.d0**dble(i)) then
              exit
           endif
        enddo
        rhomax=10.d0**dble(i)
        if (2d0*maxlargestrho.gt.rhomax) then
           rhomax=5d0*(10.d0**dble(i))
        endif
        if (verbose) then
           print *,'Automatically choosing rhomax for lammps output file:'
           print *,'Largest rho value from fitting database=',maxlargestrho
        endif
    endif
    if (verbose) then
       print *,'rhomax=',rhomax
    endif

    !Variables for potential plotting and LAMMPS potential file
    Nr=1000 !No. points used to define radial functions
    Nrho=10000 !No. points used to define embedding functions
    if ( (mod(Nr,5).ne.0).or.(mod(Nrho,5).ne.0) ) then
        print *,'ERROR:'
        print *,'   Please choose Nr and Nrho to be multiples of 5'
        print *,'to ensure correct writing of data to .eam file,'
        print *,'stopping.'
        stop
    endif
    dr=(p_rmax/dble(Nr-1)) !For LAMMPS potential file, units Angstrom
    drho=(rhomax/dble(Nrho-1)) !For LAMMPS potential file, units
    !Angstrom^-3

    !Variables for Camelion output
    !r-grid:
    !Camelion requires r-grid parameters entered into goion file in
    !format:
    !NPOT, RCUT, DELTAS, 0, 0, 0
    !Meaning: no. points on grid; radial cut-off; delta r^2
    !(first r^2 point is DELTA (rather than zero) )
    !A separate line is included for each isp, jsp, although here they
    !are defined once for all isp, jsp.
    !Note: since Camelion uses in general more radial points (it
    !doesn't use splines), we define Npot separately to Nr.
    Npot=10000 !7000
    deltar2=(p_rmax**2)/dble(Npot)
    !embfunc-grid:
    Nebf=14000!10000
    deltarho=rhomax/dble(Nebf-1) !This time first rho is zero
    delta=0.001d0 !This is use to calculate the virials

    !If MEAM potential has been optimized, write out to standard output the
    !necessary changes which need to be made to the Camelion goion file
    if ((lmax.gt.0).and.(nsteps.eq.1)) then
       print *
       print *,'Camelion output:'
       print *,'----------------'
       print *,'Substitute the following lines into the goion script (note, reduced units):'
       print *

       do isp=1,maxspecies
         do jsp=1,maxspecies
            if (jsp.ge.isp) then
               write(string2,'(A2)') element(speciestoZ(isp))
               write(string3,'(A2)') element(speciestoZ(jsp))
               string1=" set "//trim(string2)//trim(string3)//"_PV"
               !Note, in camelion, rcut does
               !not have to equal (Npot*deltar2)^0.5, it can be less. Then
               !camelion uses a simpler function for the r>rcut part.
               !      print *,Npot,p_rmax/2.684d0,deltar2/(2.684d0**2),0,0,0
               stringLen=len_trim(string1)
               if (stringLen.eq.10) then
                  write(*,'(A10,A2,I5,A1,F17.14,A1,E22.15,A13)') string1,'=(',Npot,' ', &
                      p_rmax/2.684d0,' ',deltar2/(2.684d0**2), &
                      ' 0.0 0.0 0.0)'
               elseif (stringLen.eq.11) then
                  write(*,'(A11,A2,I5,A1,F17.14,A1,E22.15,A13)') string1,'=(',Npot,' ', &
                      p_rmax/2.684d0,' ',deltar2/(2.684d0**2), &
                      ' 0.0 0.0 0.0)'
               elseif (stringLen.eq.12) then
                  write(*,'(A12,A2,I5,A1,F17.14,A1,E22.15,A13)') string1,'=(',Npot,' ', &
                      p_rmax/2.684d0,' ',deltar2/(2.684d0**2), &
                      ' 0.0 0.0 0.0)'
               endif
            endif
         enddo
       enddo
       do isp=1,maxspecies
          write(string2,'(A2)') element(speciestoZ(isp))
          string1=" set "//trim(string2)//"_EDV"
          !(the T here means 'true', in the context of angular densities)
          !      print *,'T',Npot,p_rmax/2.684d0,deltar2/(2.684d0**2),0
          stringLen=len_trim(string1)
          if (stringLen.eq.10) then
             write(*,'(A10,A4,I5,A1,F17.14,A1,E22.15,A5)') string1,'=(T ',Npot,' ', &
                 p_rmax/2.684d0,' ',deltar2/(2.684d0**2), &
                 ' 0.0)'
          elseif (stringLen.eq.11) then
             write(*,'(A11,A4,I5,A1,F17.14,A1,E22.15,A5)') string1,'=(T ',Npot,' ', &
                 p_rmax/2.684d0,' ',deltar2/(2.684d0**2), &
                 ' 0.0)'
          endif
       enddo
       do isp=1,maxspecies
          write(string2,'(A2)') element(speciestoZ(isp))
          string1=" set "//trim(string2)//"_EB"
          stringLen=len_trim(string1)
          if (stringLen.eq.9) then
             write(*,'(A9,A2,I5,A1,F17.10,A1,F17.14,A1,F17.14,A1)') string1,'=(',Nebf,' ', &
                 rhomax/0.0517d0,' ',deltarho/0.0517d0, &
                 ' ',-deltarho/0.0517d0,')'
          elseif (stringLen.eq.10) then
             write(*,'(A10,A2,I5,A1,F17.10,A1,F17.14,A1,F17.14,A1)') string1,'=(',Nebf,' ', &
                 rhomax/0.0517d0,' ',deltarho/0.0517d0, &
                 ' ',-deltarho/0.0517d0,')'
          endif
       enddo

       do isp=1,maxspecies
          write(string2,'(A2)') element(speciestoZ(isp))
          do l=1,3
             write(string3,'(I1)') l
             string1=" set "//trim(string2)//"_TWH"//trim(string3)//"="
             write(string3,'(F20.16)') meamtau(l,isp)
             string1=trim(string1)//adjustl(string3)
             write(*,'(A30)') string1
          enddo
       enddo

       print *
       write(*,*) '!!! Also change the POTENTIALDIRECTORY= line to current directory !!!'
       print *
       write(*,*) 'Please also add appropriate values for:'
       do isp=1,maxspecies
          write(string2,'(A2)') element(speciestoZ(isp))
          write(*,'(A5,A2,A31)') ' set ',string2,'_MASS_FORCECONSTANT=(#num #num)'
       enddo
       do isp=1,maxspecies
          write(string2,'(A2)') element(speciestoZ(isp))
          stringLen=len_trim(string2)
          if (stringLen.eq.1) then
             write(*,'(A11,A1,A5)') ' set OMEGA_',trim(string2),'=#num'
          elseif (stringLen.eq.2) then
             write(*,'(A11,A2,A5)') ' set OMEGA_',trim(string2),'=#num'
          endif
       enddo
       print *,'(please see Camelion documentation for further details)'

    endif
    !---------------------------------------

    if (lmax.eq.0) then
       !Write pre-amble to LAMMPS input file
       string2=""
       do i=1,maxspecies
          write(string1,'(A2)') element(speciestoZ(i))
          string2=trim(string2)//trim(string1)
       enddo
       string2=trim(string2)//" potential in LAMMPS format"
       write(58,*) trim(string2)
       write(58,*)
       write(58,*)
       if (maxspecies.lt.10) then
           write(string1,'(I1)') maxspecies
       else
           print *,'ERROR: more than 10 species in plotfunctions?!'
           print *,'STOPPING.'
           stop
       endif
       string2=trim(string1)
       do i=1,maxspecies
          write(string1,'(A2)') element(speciestoZ(i))
          string2=trim(string2)//" "//trim(string1)
       enddo
       write(58,*) trim(string2)
       write(58,'(I13,A2,F25.16,A2,I13,A2,E25.16,A2,F25.16)') Nrho,'  ',drho,'  ',Nr,'  ',dr,'  ',p_rmax
    endif

    do istr=1,nstruct
        call radialdensityfunction
        do i=1,gn_inequivalentsites(istr)
            do jj=1,gn_neighbors(i,istr)
                j=gneighborlist(jj,i,istr)
                rij=distance(i,j)
                if (distmin.eq.0d0) then
                    distmin=rij
                else
                    if (rij.lt.distmin) then
                        distmin=rij
                    endif
                endif
                if (distmax.eq.0d0) then
                    distmax=rij
                else
                    if (rij.gt.distmax) then
                        distmax=rij
                    endif
                endif
                isp=gspecies(i,istr)
                jsp=gspecies(j,istr)
                if ((isp.eq.1).and.(jsp.eq.1)) then
                    call pairpotential(isp,jsp,rij, &
                        pairpot)
                endif
            enddo
        enddo
    enddo

    gspeciesprev=gspecies(1,1)
    istr=1
    LrgVal=0d0

    do isp=1,maxspecies

       filename="embfunc"
       write(string1,'(I1)') isp
       filename=trim(filename)//trim(string1)
       open(55,file=trim(adjustl(filename)))

       filename="elecdens"
       write(string1,'(I1)') isp
       filename=trim(filename)//trim(string1)
       open(57,file=trim(adjustl(filename)))
       filename=trim(filename)//"_full"
       open(56,file=trim(adjustl(filename)))

       if (lmax.eq.0) then
          !The final two elements in the following write command are not used in
          !LAMMPS, so are placed here as place-holders.
          write(58,'(I3,A2,F15.10,A22)') speciestoZ(isp),'  ',mass(speciestoZ(isp)),' 2.85531246000000 bcc'
          !units of grams/mol (correct for 'metal' units)
       endif

       !Plot embedding function to seperate file as well as to LAMMPS file
       gspecies(1,1)=isp
       if (allocated(meam_f)) deallocate(meam_f)
       allocate(meam_f(gn_inequivalentsites(1)))
       do i=1,Nrho-4,5
           !Write out five values per line
           do j=1,5
               rho_i(1)=dble(i+j-2)*dble(rhomax)/dble(Nrho-1)
               call embeddingfunction
               meam_out(j)=meam_f(1)
               rho_out(j)=rho_i(1)
               if (all) then
                write(55,*) rho_out(j),meam_out(j)
               endif
           enddo
           if (lmax.eq.0) then
              write(58,'(e24.16,e24.16,e24.16,e24.16,e24.16)') &
                  meam_out(1),meam_out(2),meam_out(3), &
                  meam_out(4),meam_out(5)
           endif
       enddo

       !Plot electron density to seperate file as well as to LAMMPS file
       do i=1,Nr-4,5
           !Write out five values per line
           do j=1,5
               rij=dble(i+j-2)*dble(p_rmax)/dble(Nr-1)
               call raddens(rij,1,isp,density)
               density_out(j)=density(0)
               sepn_out(j)=rij
               if (density(0).ne.0d0) then
                   if (all) then
                    write(56,*) rij,density(0)
                   endif
                   if ((rij.gt.(distmin)).and. &
                       (rij.lt.(distmax))) then
                       if (all) then
                        write(57,*) rij,density(0)
                       endif
                       LrgVal=max(LrgVal,sqrt(density(0)**2))
                   endif
               endif
           enddo
           if (lmax.eq.0) then
              write(58,'(e24.16,e24.16,e24.16,e24.16,e24.16)') &
                  density_out(1),density_out(2),density_out(3), &
                  density_out(4),density_out(5)
           endif
       enddo

       close(55)
       close(56)
       close(57)

    enddo

    !Append to the end electron density files the positions of cutoffs

    do isp=1,maxspecies
       filename="elecdens"
       write(string1,'(I1)') isp
       filename=trim(filename)//trim(string1)
       open(57,file=trim(adjustl(filename)),access='append')
       do i=2*m3+lm1*m1+12*lm1*(isp-1)+2,2*m3+lm1*m1+12*lm1*(isp-1)+12,2
          if (p(i).ne.0d0) then
             if (all) then
              write(57,*)
              write(57,*) p(i),0d0
              write(57,*) p(i),LrgVal/5d0
             endif
          endif
       enddo
       close(57)
    enddo

    !Write embedding functions and electron densities to Camelion file
    if ((lmax.eq.3).and.all) then

        filename=""
        do i=1,maxspecies
            write(string1,'(A2)') element(speciestoZ(i))
            filename=trim(filename)//trim(string1)
        enddo
        string2='_MEAM.eb'
        filename=trim(filename)//trim(string2)
        open(64,file=trim(adjustl(filename)))

        do l=0,3
           filename=""
           do i=1,maxspecies
               write(string1,'(A2)') element(speciestoZ(i))
               filename=trim(filename)//trim(string1)
           enddo
           string2='_MEAM'
           filename=trim(filename)//trim(string2)
           write(string1,'(I1)') l
           filename=trim(filename)//trim(string1)
           filename=trim(filename)//".edv"
           open(59+l,file=trim(adjustl(filename)))
        enddo

        do isp=1,maxspecies
            gspecies(1,1)=isp
            istr=1
            if (allocated(meam_f)) deallocate(meam_f)
            allocate(meam_f(gn_inequivalentsites(1)))
            do i=1,Nebf
                rho_i(1)=dble(i-1)*deltarho
                call embeddingfunction
                meam_out(1)=meam_f(1)
                rho_i(1)=rho_i(1)+delta
                call embeddingfunction
                rho_i(1)=rho_i(1)-delta
                meam_out(2)=(meam_f(1)-meam_out(1))/delta
                write(64,*) meam_out(1)/0.0872d0, &
                    meam_out(2)*0.0517d0/0.0872d0
            enddo
        enddo
        gspecies(1,1)=gspeciesprev
        do l=0,3
            do isp=1,maxspecies
                do i=1,Npot
                    rij=sqrt(deltar2*dble(i))
                    call raddens(rij,1,isp,density)
                    density_out(1)=density(l)
                    rij=rij+delta
                    call raddens(rij,1,isp,density)
                    rij=rij-delta
                    density_out(2)=rij*(density(l)-density_out(1))/delta
                    !calc dens and r*d dens/dr
                    write(59+l,*) density_out(1)/0.0517d0, & !Convert to
                        density_out(2)/0.0517d0    !camelion units
                enddo

            enddo
        enddo

        close(59)
        close(60)
        close(61)
        close(62)
        close(64)

    endif

    !Plot pair-potentials seperately and to LAMMPS file
    LrgVal=0d0
    do isp=1,maxspecies
      do jsp=1,maxspecies
         if (jsp.ge.isp) then
            filename="pairpot"
            write(string1,'(I1)') isp
            filename=trim(filename)//trim(string1)
            write(string1,'(I1)') jsp
            filename=trim(filename)//trim(string1)
            open(54,file=trim(adjustl(filename)))
            filename=trim(filename)//"_full"
            open(53,file=trim(adjustl(filename)))

            do i=1,Nr-4,5
                !Write out five values per line
                do j=1,5
                    sepn=dble(i+j-2)*dble(p_rmax)/dble(Nr-1)
                   ! print *,i,'/',Nr-4,' j=',j,' sepn=',sepn,'/p',p_rmax
                    if ((sepn.eq.0d0).or.(sepn.eq.p_rmax)) then
                        pairpot=0d0 !Don't call pairpotential in this case,
                        !as the SPLINT call will fail.
                    else
                        call pairpotential(isp,jsp,sepn,pairpot)
                    endif
                    !print *,'A',pairpot
                    pairpot_out(j)=sepn*pairpot
                    sepn_out(j)=sepn
                    !print *, 'sepn=',sepn,' pairpot=',pairpot
                    if (pairpot.ne.0d0) then
                        if (all) then
                         write(53,*) sepn,pairpot
                        endif
                        !print *,'pairpot=0'
                        if ((sepn.gt.(distmin)).and. &
                            (sepn.lt.(distmax))) then
                            if (all) then
                             write(54,*) sepn,pairpot
                            endif
                            LrgVal=max(LrgVal,sqrt(pairpot**2))
                        endif
                    endif
                enddo
                if (i.eq.1) then
                    pairpot_out(1) = 2d0*pairpot_out(2) - pairpot_out(3)
                endif
                if (lmax.eq.0) then
                   write(58,'(e24.16,e24.16,e24.16,e24.16,e24.16)') &
                       pairpot_out(1),pairpot_out(2),pairpot_out(3), &
                       pairpot_out(4),pairpot_out(5)
                endif
            enddo
            close(53)
            close(54)
          endif
       enddo
    enddo

    !Append to the end of the file the positions of cutoffs
    if (all) then
       cnt=1
       do isp=1,maxspecies
          do jsp=1,maxspecies
             if (jsp.ge.isp) then
                filename="pairpot"
                write(string1,'(I1)') isp
                filename=trim(filename)//trim(string1)
                write(string1,'(I1)') jsp
                filename=trim(filename)//trim(string1)
                open(54,file=trim(adjustl(filename)),access='append')
                !Append to the end of the file the positions of cutoffs
                do i=2*m3+(4+lm1)*m1+12*lm1*m2+(cnt-1)*32+2,2*m3+(4+lm1)*m1+12*lm1*m2+cnt*32,2
                   if (p(i).ne.0d0) then
                      if (all) then
                       write(54,*)
                       write(54,*) p(i),0d0
                       write(54,*) p(i),LrgVal/5d0
                      endif
                   endif
                enddo
                close(54)
             endif
             cnt=cnt+1
          enddo
       enddo
    endif

    !Plot pair-potentials for Camelion
    if ((lmax.eq.3).and.all) then

       filename=""
       do i=1,maxspecies
          write(string1,'(A2)') element(speciestoZ(i))
          filename=trim(filename)//trim(string1)
       enddo
       string2='_MEAM.pv'
       filename=trim(filename)//trim(string2)
       open(63,file=trim(adjustl(filename)))

       do isp=1,maxspecies
          do jsp=1,maxspecies
             if (jsp.ge.isp) then

                do i=1,Npot
                   sepn=sqrt(deltar2*dble(i))
                   call pairpotential(isp,jsp,sepn,pairpot)
                   pairpot_out(1)=pairpot
                   sepn=sepn+delta
                   call pairpotential(isp,jsp,sepn,pairpot)
                   sepn=sepn-delta
                   pairpot_out(2)=sepn*(pairpot-pairpot_out(1))/delta
                   !calc pairpot and r*d pairpot/dr
                   if (all) then
                      write(63,*) pairpot_out(1)/0.0872d0, &
                          pairpot_out(2)/0.0872d0
                   endif
                enddo

             endif
          enddo
       enddo

       close(63)

    endif

    gspecies(1,1)=gspeciesprev

end subroutine plotfunctions
subroutine p_to_variables

    !---------------------------------------------------------------c
    !
    !     Copies the p() array onto the sensibly-named MEAM 
    !     variables (cmin, etc).
    !
    !     Called by:     program MEAMfit,initializemeam,meamenergy
    !     Calls:         -
    !     Arguments:     -
    !     Returns:       cmin,cmax,meamtau,
    !                 meamrhodecay,meame0,meamrho0,
    !                 meamemb3,meamemb4,pairpotparameter,
    !                 rs,rc
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff
    !
    !---------------------------------------------------------------c

    use m_filenames
    use m_meamparameters
    use m_optimization
    use m_atomproperties
    use m_generalinfo

    implicit none

    integer i,j,k,l,m

    j=1
    k=1
    l=1
    do i=1,m3
        cmin(j,k,l)=p(i)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        elseif (l.lt.maxspecies) then
            j=1
            k=1
            l=l+1
        endif
    enddo

    j=1
    k=1
    l=1
    do i=m3+1,2*m3
        cmax(j,k,l)=p(i)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        elseif (l.lt.maxspecies) then
            j=1
            k=1
            l=l+1
        endif
    enddo
    !      print *,'cmax=',cmax
    !      cmax(1:maxspecies,1:maxspecies,1:maxspecies)=p(1+m3:2*m3)

    j=1
    k=0
    do i=2*m3+1,2*m3+lm1*m1
        meamtau(k,j)=p(i)
        !print *,'Converting p to meamtau: meamtau(',k,',',j,')=',p(i),' (=p(',i,')'
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.lmax) then
            j=1
            k=k+1
        endif
    enddo
         ! print *,'meamtau=',meamtau
    !      meamtau(0:lmax,1:maxspecies)=p(2*m3+1:2*m3+lm1*m1)

    j=1
    k=0
    l=1
    m=1
    do i=2*m3+lm1*m1+1,2*m3+lm1*m1+12*lm1*m2
        meamrhodecay(j,k,l,m)=p(i)
        if (j.lt.12) then
            j=j+1
        elseif (k.lt.lmax) then
            j=1
            k=k+1
        elseif (m.lt.maxspecies) then
            j=1
            k=0
            m=m+1
        elseif (l.lt.maxspecies) then
            j=1
            k=0
            l=l+1
            m=1
        endif
    enddo
    !print *,'meamrhodecay=',meamrhodecay

    meame0(1:maxspecies)=p(2*m3+lm1*m1+12*lm1*m2+1: &
        2*m3+m1+lm1*m1+12*lm1*m2)
          !print *,'meame0=',meame0
    meamrho0(1:maxspecies)=p(2*m3+m1+lm1*m1+12*lm1*m2+1: &
        2*m3+2*m1+lm1*m1+12*lm1*m2)
          !print *,'meamrho0=',meamrho0
    meamemb3(1:maxspecies)=p(2*m3+2*m1+lm1*m1+12*lm1*m2+1: &
        2*m3+3*m1+lm1*m1+12*lm1*m2)
          !print *,'meamemb3=',meamemb3
    meamemb4(1:maxspecies)=p(2*m3+3*m1+lm1*m1+12*lm1*m2+1: &
        2*m3+4*m1+lm1*m1+12*lm1*m2)
          !print *,'meamemb4=',meamemb4

    j=1
    if (maxspecies.eq.1) then
        do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,2*m3+(4+lm1)*m1+12*lm1*m2+32
            pairpotparameter(j,1,1)=p(i)
            j=j+1
        enddo
    else
        k=1
        l=1
        do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,m4
            pairpotparameter(j,k,l)=p(i)
            if (j.lt.32) then
                j=j+1
            elseif (l.lt.maxspecies) then
                j=1
                l=l+1
            elseif (k.lt.maxspecies) then
                j=1
                l=1
                k=k+1
            endif
        enddo
    endif
    !     print *,'pairpotparameter(1:1)=',pairpotparameter(1:32,1,1)
    !     print *,'pairpotparameter(1:2)=',pairpotparameter(1:32,1,2)
    !     print *,'pairpotparameter(1:3)=',pairpotparameter(1:32,1,3)
    !     print *,'pairpotparameter(2:2)=',pairpotparameter(1:32,2,2)
    !     print *,'pairpotparameter(2:3)=',pairpotparameter(1:32,2,3)
    !     print *,'pairpotparameter(3:3)=',pairpotparameter(1:32,3,3)

    !      pairpotparameter(1:32,1,1)=p(2*m3+(4+lm1)*m1+6*lm1*m2+1:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+32)
    !      pairpotparameter(1:32,2,1)=p(2*m3+(4+lm1)*m1+6*lm1*m2+33:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+64)
    !      pairpotparameter(1:32,1,2)=p(2*m3+(4+lm1)*m1+6*lm1*m2+65:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+96)
    !      pairpotparameter(1:32,2,2)=p(2*m3+(4+lm1)*m1+6*lm1*m2+97:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+128)

    j=1
    k=1
    do i=m4+1,m4+m2
        rs(j,k)=p(i)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        endif
    enddo

    j=1
    k=1
    do i=m4+m2+1,m4+2*m2
        rc(j,k)=p(i)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        endif
    enddo

    j=1
    do i=m4+2*m2+1,m4+2*m2+m1
        enconst(j)=p(i)
        j=j+1
    enddo
    !print *,'enconst=',enconst

end subroutine p_to_variables
function pythag(a,b)
    implicit none
    real(8) a,b,pythag,absa,absb
    absa=abs(a)
    absb=abs(b)
    if(absa.gt.absb)then
        pythag=absa*sqrt(1d0+(absb/absa)**2)
    else
        if(absb.eq.0d0)then
            pythag=0d0
        else
            pythag=absb*sqrt(1d0+(absa/absb)**2)
        endif
    endif
end function pythag
subroutine raddens(rij,i,j,density)

    use m_electrondensity
    use m_meamparameters

    implicit none

    integer i,j,n,o
    real(8) rij,density(0:lmax)

    !     !---- Comment in the following block for Hepburn potential ----
    !     real(8) r1CFel0,r2CFel0,r3,r4,
    !    +     r1FeCl0,r2FeCl0,
    !    +     A_CFe1l0,B_CFe1l0,C_CFe1l0,D_CFe1l0,
    !    +     A_CFe2,B_CFe2,C_CFe2,D_CFe2,
    !    +     A_FeCl0,B_FeCl0,C_FeCl0,D_FeCl0,
    !    +     denom,thi_sat_l0,tmp
    !
    !     if (((i.eq.1).and.(j.eq.1)).or.
    !    +     (i.eq.2).and.(j.eq.2)) then
    !        density(0) = 0d0
    !        if (rij.le.meamrhodecay(2,0,i,j)) then
    !            density(0) = density(0) +
    !    +           meamrhodecay(1,0,i,j) *
    !    +           ((meamrhodecay(2,0,i,j) - rij) **3)
    !        endif
    !        if (rij.le.meamrhodecay(4,0,i,j)) then
    !            density(0) = density(0) +
    !    +            meamrhodecay(3,0,i,j) *
    !    +           ((meamrhodecay(4,0,i,j) - rij) **3)
    !        endif
    !        if (rij.le.meamrhodecay(6,0,i,j)) then
    !            density(0) = density(0) +
    !    +            meamrhodecay(5,0,i,j) *
    !    +            ((meamrhodecay(6,0,i,j) - rij) **3)
    !        endif
    !        if (rij.le.meamrhodecay(8,0,i,j)) then
    !            density(0) = density(0) +
    !    +            meamrhodecay(7,0,i,j) *
    !    +            ((meamrhodecay(8,0,i,j) - rij) **3)
    !        endif
    !        if (rij.le.meamrhodecay(10,0,i,j)) then
    !            density(0) = density(0) +
    !    +            meamrhodecay(9,0,i,j) *
    !    +            ((meamrhodecay(10,0,i,j) - rij) **3)
    !        endif
    !        if (rij.le.meamrhodecay(12,0,i,j)) then
    !            density(0) = density(0) +
    !    +            meamrhodecay(11,0,i,j) *
    !    +           ((meamrhodecay(12,0,i,j) - rij) **3)
    !        endif
    !     elseif ((i.eq.2).and.(j.eq.1)) then
    !        r1CFel0=min(max(meamrhodecay(1,0,2,1),1.70d0),1.85d0)
    !        r2CFel0=min(max(meamrhodecay(2,0,2,1),1.71d0),1.86d0)
    !        r3=         meamrhodecay(3,0,2,1)
    !        r4=         meamrhodecay(4,0,2,1)
    !        denom=3d0*(r1CFel0+r2CFel0)*
    !    +       ((r1CFel0**2)-(r2CFel0**2))
    !    +       - 4d0*((r1CFel0**3)-(r2CFel0**3))
    !        A_CFe1l0=0.001d0 + (0.5d0-0.001d0) *
    !    +      (-3d0*r1CFel0*(r2CFel0**2) + (r2CFel0**3))
    !    +      / denom
    !        B_CFe1l0=6d0*(0.5d0-0.001d0)*r1CFel0*r2CFel0
    !    +      / denom
    !        C_CFe1l0=-3d0*(0.5d0-0.001d0)*
    !    +      (r1CFel0+r2CFel0) / denom
    !        D_CFe1l0=2d0*(0.5d0-0.001d0) / denom
    !        denom=3d0*(r3+r4)*((r3**2)-(r4**2))
    !    +       - 4d0*((r3**3)-(r4**3))
    !        A_CFe2=0.001d0 *
    !    +      (-3d0*r3*(r4**2) + (r4**3)) / denom
    !        B_CFe2=6d0*0.001d0*r3*r4 / denom
    !        C_CFe2=-3d0*0.001d0*(r3+r4) / denom
    !        D_CFe2=2d0*0.001d0 / denom
    !        if (rij.le.r1CFel0) then
    !           density(0) = 0.5d0
    !        elseif ((rij.gt.r1CFel0).and.
    !    +          (rij.le.r2CFel0)) then
    !           density(0) = A_CFe1l0 + B_CFe1l0*rij
    !    +                    + C_CFe1l0*(rij**2) + D_CFe1l0*(rij**3)
    !        elseif ((rij.gt.r2CFel0).and.(rij.le.r3)) then
    !           density(0) = 0.001d0
    !        elseif ((rij.gt.r3).and.(rij.le.r4)) then
    !           density(0) = A_CFe2 + B_CFe2*rij
    !    +                    + C_CFe2*(rij**2) + D_CFe2*(rij**3)
    !        elseif (rij.gt.r4) then
    !           density(0) = 0d0
    !        endif
    !     elseif ((i.eq.1).and.(j.eq.2)) then
    !        r1FeCl0=meamrhodecay(1,0,1,2)
    !        r2FeCl0=meamrhodecay(2,0,1,2)
    !        thi_sat_l0= meamrhodecay(3,0,1,2)
    !        denom=3d0*(r1FeCl0+r2FeCl0)*((r1FeCl0**2)-
    !    +       (r2FeCl0**2))- 4d0*
    !    +       ((r1FeCl0**3)-(r2FeCl0**3))
    !        A_FeCl0= thi_sat_l0 *
    !    +      (-3d0*r1FeCl0*(r2FeCl0**2) + (r2FeCl0**3))
    !    +      / denom
    !        B_FeCl0=6d0*thi_sat_l0*r1FeCl0*r2FeCl0 / denom
    !        C_FeCl0=-3d0*thi_sat_l0*(r1FeCl0+r2FeCl0)
    !    +      / denom
    !        D_FeCl0=2d0*thi_sat_l0 / denom
    !        if (rij.le.r1FeCl0) then
    !           density(0) = thi_sat_l0
    !        elseif ((rij.gt.r1FeCl0).and.
    !    +           (rij.le.r2FeCl0)) then
    !           density(0) = A_FeCl0 + B_FeCl0*rij
    !    +                    + C_FeCl0*(rij**2) + D_FeCl0*(rij**3)
    !        elseif (rij.gt.r2FeCl0) then
    !           density(0) = 0d0
    !        endif
    !     endif
    !     density(1)=0d0
    !     density(2)=0d0
    !     density(3)=0d0
    !     !--------------------------------------------------------------

    !---- Comment in the following block (and comment ----
    !     out the other blocks), to use the sum over
    !     cubic-term form of the radial density
    do o=0,lmax !loop over ang. mom.
        density(o) = 0d0
        do n=2,12,2
            if (rij.le.meamrhodecay(n,o,i,j)) then
                density(o) = density(o) + &
                    meamrhodecay(n-1,o,i,j) * &
                    ((meamrhodecay(n,o,i,j) - rij) **3)
            endif
        enddo
    enddo
    !      print *,'density(1)=',density(1)
    !      print *,'meamrhodecay(1:4,1,',i,',',j,')=',
    !     +         meamrhodecay(1:4,1,i,j)
    !      print *,'rij=',rij
    !      stop

    !      if ((density(0).gt.0.523918344152374d0).and.
    !     +       (density(0).lt.0.523918344152375d0)) then
    !         write(*,*) density(0),i,j,rij
    !         print *,meamrhodecay(1:12,0,i,j)
    !         stop
    !      endif


    !      print *,'meamrhodecay(1:12,0,',i,',',j,')=',
    !     +         meamrhodecay(1:12,0,i,j)
    !      print *,density
    !c      stop
    !-----------------------------------------------------


    !      !---- Comment in the following block (and comment ----
    !      !     out the other blocks), to use the (hard-
    !      !     wired) Lau 2007 potential
    !      if ((i.eq.1).and.(j.eq.1)) then
    !         !Fe-Fe
    !         if ((rij.le.3.569745d0)) then
    !            density(0) = ((rij-3.569745d0)**2)
    !     +         + 0.504238d0 * ((rij-3.569745d0)**3)
    !         else
    !            density(0) = 0d0
    !         endif
    !      elseif ((i.eq.2).and.(j.eq.1)) then
    !         if (rij.lt.2.545937d0) then
    !            density(0) =
    !     +         10.482408d0 * ((rij-2.545937d0)**2)
    !     +         + 3.782595d0 * ((rij-2.545937d0)**3)
    !         else
    !            density(0) = 0d0
    !         endif
    !      elseif ((i.eq.1).and.(j.eq.2)) then
    !         if (rij.lt.2.545937d0) then
    !            density(0) =
    !     +         10.024001d0 * ((rij-2.545937d0)**2)
    !     +         + 1.638980d0 * ((rij-2.545937d0)**3)
    !         else
    !            density(0) = 0d0
    !         endif
    !      elseif ((i.eq.2).and.(j.eq.2)) then
    !         if (rij.lt.2.892070d0) then
    !            density(0) =
    !     +         - 7.329211d0 * ((rij-2.892070d0)**3)
    !         else
    !            density(0) = 0d0
    !         endif
    !      endif
    !      density(1) = 0d0
    !      density(2) = 0d0
    !      density(3) = 0d0
    !      !-----------------------------------------------------

end subroutine raddens
subroutine radialcutofffunction(q,hc)

    !------------------------------------------------------------c
    !
    !     Calculates the radial cutoff function, hc, used for
    !     radial density function. hc=1 if q<0; hc decreases
    !     from 1 to 0 over the interval 0<q<1; hc=0 if q>1.
    !
    !     Called by:     radialdensityfunction
    !     Calls:         -
    !     Arguments:     q
    !     Returns:       hc
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter, Feb 9 2006
    !
    !------------------------------------------------------------c

    implicit none

    real*8 q,hc

    hc=1d0
    if(q.ge.1d0) then
        hc=0d0
    else if(q.gt.0d0) then
        hc=1d0-q*q*q*(6d0*q*q-15d0*q+10d0)
    endif

end subroutine radialcutofffunction
subroutine radialdensityfunction

    !----------------------------------------------------------------c
    !
    !     Fills the module m_radialdensity with the f_j_L_ij
    !     parameters. NOTE: I use f_j_L_ij = rhoatnucleus_j *
    !     exp( -decay(i,j,L)*r_ij )
    !
    !     Called by:     meamenergy
    !     Calls:         distance,radialcutofffunction
    !     Arguments:     istr,gn_inequivalentsites,gn_neighbors,
    !                 lmax,gspecies,gn_neighbors,gneighborlist,rs,
    !                 rc,meamrhodecay
    !     Returns:       fjlij
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter, Feb 9 2006
    !
    !----------------------------------------------------------------c

    use m_generalinfo
    use m_atomproperties
    use m_geometry
    use m_meamparameters
    use m_electrondensity   !f_j_L_ij

    implicit none

    integer i,j,jj,nni,iaux(1),isp,jsp,n,o

    if(allocated(fjlij)) deallocate(fjlij)
    iaux=maxval(gn_neighbors(1:gn_inequivalentsites(istr),istr))
    allocate(fjlij(0:lmax,iaux(1),gn_inequivalentsites(istr)))

    do i=1,gn_inequivalentsites(istr)
        isp=1
        if (thiaccptindepndt.eqv..false.) isp=gspecies(i,istr)
        nni=gn_neighbors(i,istr)

        do jj=1,nni
            j=gneighborlist(jj,i,istr)
            jsp=gspecies(j,istr)
            rij=diststr(jj,i,0,0,istr)
            !            rij=distance(i,j)
            !           if (rij.ne.diststr(jj,i,0,0,istr)) then
            !              print *,'rij does not match the stored value'
            !              print *,'istr=',istr,', i=',i,',j=',j
            !              print *,'rij calc=',rij,', stored=',
            !     +    diststr(jj,i,0,0,istr)
            !              stop
            !           endif
            call raddens(rij,isp,jsp,fjlij(0:lmax,jj,i))
        enddo

    enddo

   !if (istr.eq.1) then
   !do j=1,m1
   !do o=0,lmax !loop over ang. mom.
   !    do n=1,12
   !        write(81,'(A12,I2,A3,I2,A3,I2,A3,I2,A2,F15.10)') 'meamrodecay(',n,',l=',o,',i=',1,',j=',j,')=',meamrhodecay(n,o,1,j)
   !    enddo
   !enddo
   !enddo
   !endif

end subroutine radialdensityfunction
subroutine radialdensityfunction_onlyiatom(iatom,cart)

    use m_generalinfo
    use m_atomproperties
    use m_geometry
    use m_meamparameters
    use m_electrondensity   !f_j_L_ij

    implicit none

    integer i,j,jj,nni,iaux(1),isp,jsp,iatom,cart

    iaux=maxval(gn_neighbors(1:gn_inequivalentsites(istr),istr))

    isp=gspecies(iatom,istr)
    nni=gn_neighbors(iatom,istr)

    if (fastforce) then

       do jj=1,nni
           j=gneighborlist(jj,iatom,istr)
           jsp=gspecies(j,istr)
           rij=diststr(jj,iatom,iatom,cart,istr)
           !         if (rij.ne.diststr(jj,iatom,iatom,cart,istr)) then
           !            print *,'rij does not match the stored value'
           !            print *,'istr=',istr,', iatom=',iatom,',j=',j
           !            print *,'rij calc=',rij,
           !     +        ', stored=',diststr(jj,iatom,iatom,cart,istr)
           !            stop
           !         endif
    
           if (thiaccptindepndt.eqv..true.) then
               call raddens(rij,1,jsp,fjlij(0:lmax,jj,iatom))
           else
               call raddens(rij,isp,jsp,fjlij(0:lmax,jj,iatom))
           endif
       enddo
    
       do i=1,gn_inequivalentsites(istr)
    
           if (i.ne.iatom) then
               isp=gspecies(i,istr)
               nni=gn_neighbors(i,istr)
    
               do jj=1,nni
                   j=gneighborlist(jj,i,istr)
                   if (j.eq.iatom) then
                       jsp=gspecies(j,istr)
                       rij=diststr(jj,i,iatom,cart,istr)
                       !                  if (rij.ne.diststr(jj,i,iatom,cart,istr)) then
                       !                     print *,'rij does not match the stored value'
                       !                     print *,'jj=',jj,', i=',i,', istr=',istr,
                       !     +       ', iatom=',iatom,',cart=',cart
                       !                     print *,'rij calc=',rij,
                       !     +                 ', stored=',diststr(jj,i,iatom,cart,istr)
                       !                     stop
                       !                  endif
    
                       if (thiaccptindepndt.eqv..true.) then
                           call raddens(rij,1,jsp,fjlij(0:lmax,jj,i))
                       else
                           call raddens(rij,isp,jsp,fjlij(0:lmax,jj,i))
                       endif
                   endif
               enddo
    
           endif
    
       enddo

    else

       do jj=1,nni
           j=gneighborlist(jj,iatom,istr)
           jsp=gspecies(j,istr)
           rij=distance(iatom,j)
    
           if (thiaccptindepndt.eqv..true.) then
               call raddens(rij,1,jsp,fjlij(0:lmax,jj,iatom))
           else
               call raddens(rij,isp,jsp,fjlij(0:lmax,jj,iatom))
           endif
       enddo
    
       do i=1,gn_inequivalentsites(istr)
    
           if (i.ne.iatom) then
               isp=gspecies(i,istr)
               nni=gn_neighbors(i,istr)
    
               do jj=1,nni
                   j=gneighborlist(jj,i,istr)
                   if (j.eq.iatom) then
                       jsp=gspecies(j,istr)
                       rij=distance(i,j)
    
                       if (thiaccptindepndt.eqv..true.) then
                           call raddens(rij,1,jsp,fjlij(0:lmax,jj,i))
                       else
                           call raddens(rij,isp,jsp,fjlij(0:lmax,jj,i))
                       endif
                   endif
               enddo
    
           endif
    
       enddo

    endif

end subroutine radialdensityfunction_onlyiatom
subroutine readAtomPropsFromVasprun

    !----------------------------------------------------------------------c
    !
    !     From the start position <modeling> in vasprun (read in through
    !     input channel 1, extract the following: nspecies, natoms, z, 
    !     nat, zz, and tr.
    !
    !     Called by:     
    !     Calls:         
    !     Arguments:    
    !     Returns:    
    !     Files read:   
    !     Files written:
    !
    !     Andrew Duff, 2014
    !
    !----------------------------------------------------------------------c

    use m_poscar

    implicit none

    integer wordLength,i,j,ispc,it1,marker
    character*2 ChemSym
    character*80 strng1,strng2
    character*80 line
    character*80, allocatable:: words(:)

    !natoms is found under a line such as:   <atoms>      64</atoms>
    do
      read(1,'(a80)') line
      if (line(1:9).eq.'  <atoms>') then
         read(line,*) strng1,strng2
         call keepNumber(strng2)
         strng2=trim(strng2)
         !wordLength=len_trim(strng2)
         !strng2=strng2(1:wordLength-1)
         read(strng2,*) natoms
         !print *,'natoms=',natoms
         exit
      endif
    enddo

    !nspecies is found under a line such as:  <types>       2</types>
    do
      read(1,'(a80)') line
      if (line(1:9).eq.'  <types>') then
         read(line,*) strng1,strng2
         call keepNumber(strng2)
         strng2=trim(strng2)
         !wordLength=len_trim(strng2)
         !strng2=strng2(1:wordLength-1)
         read(strng2,*) nspecies
         !print *,'nspecies=',nspecies
         exit
      endif
    enddo
 
    if (.not.allocated(nat)) then
       allocate(nat(nspecies),z(nspecies))
       allocate(zz(natoms))
    endif
 
    do i=1,5
       read(1,*)
    enddo
    !atom types are of the form:     <rc><c>Zr</c><c>   1</c></rc>
 
    allocate(words(nspecies))
 
    read(1,*) strng1
    words(1)=strng1(8:9)
    !print *,'first chem symbol=',words(1)
    ispc=1
    marker=1
    do i=2,natoms
       read(1,*) strng1
       !print *,'read strng1=',strng1
       ChemSym(1:2)=strng1(8:9)
       !print *,'... atom ',i,' chem symbol=',ChemSym(1:2)
       if (ChemSym.ne.words(ispc)) then
          nat(ispc)=i-marker
          marker=i
          ispc=ispc+1
          words(ispc)=ChemSym
       endif
    enddo
    nat(ispc)=natoms-marker+1
    !print *,'nat(1:',nspecies,')=',nat(1:nspecies)
 
    !print *,'chemsymbs:',words(1:nspecies)
    call elementToZ(words)
    !print *,'Z(1:',nspecies,')=',Z(1:nspecies)
    !print *,'species per atom (zz):'
    it1=0
    do i=1,nspecies
        do j=1,nat(i)
            it1=it1+1
            zz(it1)=z(i)
            !print *,'zz(',it1,')=',zz(it1),'=z(',i,')=',z(i)
        enddo
    enddo
    !print *,'size(zz)=',Size(zz)
    !stop
 
    !Lattice vectors at: "   <varray name="basis" >"
    !print *,'trying to find <varray name, now:'
    do
      read(1,'(a80)') line
      if ((line(1:16).eq.'   <varray name=').and.(line(18:22).eq.'basis')) then
         !print *,'found <varray name= basis >'
         read(1,*) strng1,tr(1,1),tr(2,1),strng2
         !print *,'read line 1'
         call keepNumber(strng2)
         strng2=trim(strng2)
         !wordLength=len_trim(strng2)
         !strng2=strng2(1:wordLength-1)
         read(strng2,*) tr(3,1)
         read(1,*) strng1,tr(1,2),tr(2,2),strng2
         !print *,'read line 2'
         call keepNumber(strng2)
         strng2=trim(strng2)
         !wordLength=len_trim(strng2)
         !strng2=strng2(1:wordLength-1)
         read(strng2,*) tr(3,2)
         read(1,*) strng1,tr(1,3),tr(2,3),strng2
         !print *,'read line 3'
         call keepNumber(strng2)
         strng2=trim(strng2)
         !wordLength=len_trim(strng2)
         !strng2=strng2(1:wordLength-1)
         read(strng2,*) tr(3,3)
         !print *,'read line 4'
         exit
      endif
    enddo
    !print *,'finished read sub'
end subroutine readAtomPropsFromVasprun

subroutine keepNumber(strng)

!   Keeps numbers, d's, e's, g's and '-' signs in the char*80 strng. 
!   For keeping numerical components of strings (d for double precision, etc)

    implicit none

    integer i
    character*80 strng

    !print *,'starting string=',strng

    do i=1,80
       if ((strng(i:i).eq.'0').or. &
           (strng(i:i).eq.'1').or. &
           (strng(i:i).eq.'2').or. &
           (strng(i:i).eq.'3').or. &
           (strng(i:i).eq.'4').or. &
           (strng(i:i).eq.'5').or. &
           (strng(i:i).eq.'6').or. &
           (strng(i:i).eq.'7').or. &
           (strng(i:i).eq.'8').or. &
           (strng(i:i).eq.'9').or. &
           (strng(i:i).eq.'-').or. &
           (strng(i:i).eq.'d').or. &
           (strng(i:i).eq.'D').or. &
           (strng(i:i).eq.'e').or. &
           (strng(i:i).eq.'E').or. &
           (strng(i:i).eq.'g').or. &
           (strng(i:i).eq.'G').or. &
           (strng(i:i).eq.'.')) then
          !OK
       else
          strng(i:i)=' '
       endif
    enddo

    !print *,'finishing string=',strng
    !stop

end subroutine KeepNumber
subroutine readBoundsForRandParas

    !--------------------------------------------------------------c
    !
    !     Reads in the limits which are to be used to initialize
    !     the random starting values of the MEAM parameters from
    !     the 'dataforMEAMparagen' file.
    !
    !     Called by:
    !     Calls:
    !     Arguments: maxspecies,lmax
    !     Returns:
    !     Files read:
    !     Files written:
    !
    !     Andy Duff
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_filenames
    use m_atomproperties
    use m_generalinfo
    use m_optimization

    implicit none

    integer i,j,ii,ll,ispc,jspc,spc_cnt
    character*80 string
    !Variables just used for reading in start_meam_parameters file:
    integer nfreep

    allocate( meamtau_minorder(0:lmax,maxspecies), &
        meamtau_maxorder(0:lmax,maxspecies), &
        meamrhodecay_minradius(0:lmax,maxspecies,maxspecies), &
        meamrhodecay_maxradius(0:lmax,maxspecies,maxspecies), &
        meamrhodecay_negvals(0:lmax,maxspecies,maxspecies), &
        meamrhodecay_minorder(1:6,0:lmax,maxspecies,maxspecies), &
        meamrhodecay_maxorder(1:6,0:lmax,maxspecies,maxspecies), &
        meame0_minorder(maxspecies), meame0_maxorder(maxspecies), &
        meamrho0_minorder(maxspecies), meamrho0_maxorder(maxspecies), &
        meamemb3_minorder(maxspecies), meamemb3_maxorder(maxspecies), &
        meamemb4_minorder(maxspecies), meamemb4_maxorder(maxspecies), &
        pairpotparameter_minradius(maxspecies,maxspecies), &
        pairpotparameter_maxradius(maxspecies,maxspecies), &
        pairpotparameter_negvals(maxspecies,maxspecies), &
        pairpotparameter_minorder(1:8,maxspecies,maxspecies), &
        pairpotparameter_maxorder(1:8,maxspecies,maxspecies), &
        minordervaluepairpot(1:16), &
        maxordervaluepairpot(1:16), &
        nfreeppairpot(maxspecies,maxspecies), &
        enconst_minorder(maxspecies), enconst_maxorder(maxspecies) )

    !Set the maxradius variables to zero (then only those 'in use' will be used
    !in the 'checkParas' subroutine).
    meamrhodecay_maxradius=0d0
    pairpotparameter_maxradius=0d0

    open(unit=1,file=trim(settingsfile))

    do
       read(1,*) string
       if (string.eq."RANDPARASMANUAL") then
          exit
       endif
    enddo

    !meamtau
    do ll=1,lmax
        do ispc=1,maxspecies
            if (freep(2*m3+ispc+ll*m1).ge.2) then
                read(1,*) meamtau_minorder(ll,ispc)
                read(1,*) meamtau_maxorder(ll,ispc)
            endif
        enddo
    enddo

    !electron density functions
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    !Check if we need to generate random paramaters for
                    !thi(ispc,jspc,ll)
                    nfreep=0
                    do i=1,12
                        if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i).ge.2) then
                            nfreep=nfreep+1
                        endif
                    enddo
                    if (nfreep.gt.0) then
                       print *,'rho(',ispc,',',jspc,'):'
                        if (typethi(ispc,jspc).eq.1) then !Cubic form
                            !Read in parameters from 'dataforMEAMparagen'
                            !necessary
                            !to generate these parameters
                            read(1,*) meamrhodecay_maxradius(ll,ispc,jspc)
                            read(1,*) meamrhodecay_minradius(ll,ispc,jspc)
                            read(1,*) meamrhodecay_negvals(ll,ispc,jspc)
                            print *,'rho maxrad=',meamrhodecay_maxradius(ll,ispc,jspc)
                            print *,'rho minrad=',meamrhodecay_minradius(ll,ispc,jspc)
                            print *,'rho negvasl=',meamrhodecay_negvals(ll,ispc,jspc)

                            !For each coefficient which needs to be generated,
                            !read in
                            !min and max orders (10^..)
                            ii=1
                            do i=1,6
                                if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+2*i-1).ge.2) &
                                    then
                                    read(1,*) meamrhodecay_minorder(ii,ll,ispc,jspc)
                                    read(1,*) meamrhodecay_maxorder(ii,ll,ispc,jspc)
                                    print *,'rho minorder=',meamrhodecay_minorder(ii,ll,ispc,jspc)
                                    print *,'rho maxorder=',meamrhodecay_maxorder(ii,ll,ispc,jspc)
                                    ii=ii+1
                                endif
                            enddo
                        endif
                    endif
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !embedding functions
    read(1,*) embfuncRandType
    do ispc=1,maxspecies
        do i=1,4
            if (freep(2*m3+lm1*m1+12*lm1*m2+ispc+(i-1)*m1).ge.2) then
                if (i.eq.1) then
                    read(1,*) meame0_minorder(ispc)
                    read(1,*) meame0_maxorder(ispc)
                   print *,'meame0 minO=',meame0_minorder(ispc)
                   print *,'meame0 maxO=',meame0_maxorder(ispc)
                elseif (i.eq.2) then
                    read(1,*) meamrho0_minorder(ispc)
                    read(1,*) meamrho0_maxorder(ispc)
                   print *,'meamrho0 minO=',meamrho0_minorder(ispc)
                   print *,'meamrho0 maxO=',meamrho0_maxorder(ispc)
                elseif (i.eq.3) then
                    read(1,*) meamemb3_minorder(ispc)
                    read(1,*) meamemb3_maxorder(ispc)
                   print *,'meamemb3 minO=',meamemb3_minorder(ispc)
                   print *,'meamemb3 maxO=',meamemb3_maxorder(ispc)
                elseif (i.eq.4) then
                    read(1,*) meamemb4_minorder(ispc)
                    read(1,*) meamemb4_maxorder(ispc)
                   print *,'meamemb4 minO=',meamemb4_minorder(ispc)
                   print *,'meamemb4 maxO=',meamemb4_maxorder(ispc)
                endif
            endif
        enddo
    enddo

    !pair-potentials
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                !Check if we need to generate random paramaters for
                !pairpot(isp,jspc)
                nfreeppairpot(ispc,jspc)=0
                do i=1,32
                    if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt) &
                        .ge.2) then
                        nfreeppairpot(ispc,jspc)=nfreeppairpot(ispc,jspc)+1
                    endif
                enddo
                if (nfreeppairpot(ispc,jspc).gt.0) then
                    if (typepairpot(ispc,jspc).eq.2) then
                        !Read in parameters from 'dataforMEAMparagen' necessary
                        !to
                        !generate these parameters
                        read(1,*) pairpotparameter_maxradius(ispc,jspc)
                        read(1,*) pairpotparameter_minradius(ispc,jspc)
                        read(1,*) pairpotparameter_negvals(ispc,jspc)
                        print *,'pairpot maxrad=',pairpotparameter_maxradius(ispc,jspc)
                        print *,'pairpot minrad=',pairpotparameter_minradius(ispc,jspc)

                        if ((pairpotparameter_negvals(ispc,jspc).lt.1).or. &
                            (pairpotparameter_negvals(ispc,jspc).gt.3)) then
                            print *,'      ERROR: you must specify 1, 2 or 3,'
                            print *,'      stopping.'
                            stop
                        endif
                        !For each coefficient which needs to be generated, read
                        !in min and max orders (10^..)
                        ii=1
                        do i=1,16
                            if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+2*i-1+32*spc_cnt) &
                                .ge.2) then
                                read(1,*) pairpotparameter_minorder(ii,ispc,jspc)
                                read(1,*) pairpotparameter_maxorder(ii,ispc,jspc)
                                ii=ii+1
                            endif
                        enddo
                    else
                        print *,'typepairpot=',typepairpot(ispc,jspc), &
                            'not supported,stopping'
                        stop
                    endif
                endif
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !enconst
    do ispc=1,maxspecies
        if (freep(m4+2*m2+ispc).ge.2) then
            read(1,*) enconst_minorder(ispc)
            read(1,*) enconst_maxorder(ispc) 
        endif
    enddo

    close(1)

end subroutine readBoundsForRandParas
subroutine readdatapoints

    !----------------------------------------------------------------------c
    !
    !     Reads in energies and/or forces from outcar files and stores them
    !     in the 'truedata' array.
    !
    !     Called by:     program MEAMfit
    !     Calls:         readoutcarforces,readoutcarenergy
    !     Returns:       ndatapoints,truedata
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Jan 2008
    !
    !----------------------------------------------------------------------c

    use m_filenames
    use m_datapoints
    use m_optimization
    use m_geometry

    implicit none

    character*80 outcarfile_prev
    character*11 string
    integer istruc,idatapoint,nconfig_prev,optimizeforce_prev
    real(8), allocatable:: truedata_readin(:)

    print *,'Preparing to read data from vasprun.xml files'

    !Determine the number of datapoints to read in from the outcar
    !files. This will depend on whether we are fitting to energies or
    !forces for each of the outcar files.
    ndatapoints=0
    do istruc=1,nstruct
        if (optimizeforce(istruc).gt.0) then
            ndatapoints=ndatapoints+3*gn_forces(istruc)
        else
            ndatapoints=ndatapoints+1
        endif
    enddo
    !Allocate fitting arrays.
    if (allocated(truedata)) deallocate(truedata)
    if (allocated(fitdata)) deallocate(fitdata)
    if (allocated(sumppStr)) deallocate(sumppStr)
    if (allocated(summeamfStr)) deallocate(summeamfStr)
    if (allocated(bestfitdata)) deallocate(bestfitdata)
    if (allocated(smallestsepnStr)) deallocate(smallestsepnStr)
    if (allocated(smlsepnperspcStr)) deallocate(smlsepnperspcStr)
    if (allocated(lrgsepnperspcStr)) deallocate(lrgsepnperspcStr)

    allocate(truedata(ndatapoints),fitdata(ndatapoints), &
        bestfitdata(ndatapoints),sumppStr(ndatapoints), &
        summeamfStr(ndatapoints),smallestsepnStr(ndatapoints), &
        smlsepnperspcStr(10,ndatapoints), &
        lrgsepnperspcStr(10,ndatapoints) )

    !Fill out the 'truedata' array by reading the true energies/
    !forces from the vasprun.xml/OUTCAR files.
    idatapoint=1
    optimizeforce_prev=0
    do istruc=1,nstruct !Loop through all structures present in the
                        !vasprun.xml/OUTCAR files.

        if (optimizeforce(istruc).gt.0) then
           !Fit to forces:
           allocate(truedata_readin(1:3*gn_forces(istruc)))
           !First check if forces to be read in follow on immediately from the
           !previous set of forces read in (avoids re-reading entire vasprun.xml
           !file)
           if ((optimizeforce(istruc).eq.optimizeforce_prev).and. &
               (trim(outcarfiles(istruc)).eq.outcarfile_prev).and. &
               (nconfig(istruc).gt.nconfig_prev).and. &
               (istruc.gt.1)) then
              call readoutcarforces(outcarfiles(istruc), &
                  gn_forces(istruc),truedata_readin,istruc,nconfig(istruc), &
                  vasprun(istruc),nconfig_prev+1)
           else
               call readoutcarforces(outcarfiles(istruc), &
                   gn_forces(istruc),truedata_readin,istruc,nconfig(istruc), &
                   vasprun(istruc),0)
           endif
            truedata(idatapoint:idatapoint+ &
                3*gn_forces(istruc)-1)= &
                truedata_readin(1:3*gn_forces(istruc))
            deallocate(truedata_readin)
            idatapoint=idatapoint+3*gn_forces(istruc)
            !            print *,'ndatapoints=',idatapoint,' (force)'
            !            print *,'   (gn_forces(',i,')=',gn_forces(i)
        else
          if ((optimizeforce(istruc).eq.optimizeforce_prev).and. &
              (trim(outcarfiles(istruc)).eq.outcarfile_prev).and. &
              (nconfig(istruc).gt.nconfig_prev).and. &
              (istruc.gt.1)) then
              !print *,'calling for continuation readin, since'
              !print *,'nconfig(',istruc,')=',nconfig(istruc),'='
              !print *,'nconfig_prev+1'
              call readoutcarenergy( trim(outcarfiles(istruc)), &
                  truedata(idatapoint),freeenergy(istruc),nconfig(istruc), &
                  vasprun(istruc),nconfig_prev+1)
          else
              !print *,'no continuation readin, since'
              !print *,'nconfig(',istruc,')=',nconfig(istruc),'='
              !print *,'nconfig_prev+1'
               call readoutcarenergy( trim(outcarfiles(istruc)), &
                   truedata(idatapoint),freeenergy(istruc),nconfig(istruc), &
                   vasprun(istruc),0)
          endif
          !print *,'ndatapoints=',idatapoint,' energy=',truedata(idatapoint)
          idatapoint=idatapoint+1

        endif
        optimizeforce_prev=optimizeforce(istruc)
        outcarfile_prev=trim(outcarfiles(istruc))
        nconfig_prev=nconfig(istruc)

    enddo
    close(1)

    print *,'Completed data read-in from vasprun.xml files'

end subroutine readdatapoints
        subroutine readFitdbseLine(line,string1,string2,string3,string4,double1)

        implicit none

        logical spcPrev
        integer iprev,iword,i
        real(8) double1
        character*80 line
        character*80 str,string1,string2,string3,string4

        spcPrev=.true. !If this is true, don't accept the space as delimiter
        iprev=1 !Marks the start of a string to be extracted
        iword=1

        do i=1,80
          if (line(i:i).eq." ") then

            if (spcPrev.eqv..false.) then
               if (iword.eq.1) then
                  string1=trim(line(iprev:i))
               !To reinstate POSCAR/OUTCAR pairs, re-comment in
               !following lines:
!              elseif (iword.eq.2) then
!                 string2=trim(line(iprev:i))
               elseif (iword.eq.2) then
                  string3=trim(line(iprev:i))
               elseif (iword.eq.3) then
                  string4=trim(line(iprev:i))
               elseif (iword.eq.4) then
                  str=trim(line(iprev:i))
                  read( str, * ) double1
               endif
               spcPrev=.true.
               iword=iword+1
               if (iword.eq.5) then
                  exit
               endif
            endif
            spcPrev=.true.

          else
             if (spcPrev.eqv..true.) then
                iprev=i
             endif
             spcPrev=.false.
          endif
        enddo
        if (iword.ne.5) then
           print *,'ERROR: expected 4 words in fitdbse file, found only',iword-1
           if (i.eq.81) then
              print *,'(you used at least 80 characters in the line I just read...'
              print *,'please make sure to use NO MORE than 80)'
           endif
           print *,'STOPPING.'
           stop
        endif

        end subroutine readFitdbseLine

subroutine readinoptparasonebyone

    !-------------------------------------------------------------------c
    !
    !     Read in the which potential parameters are to be optimized.
    !
    !     Called by:     readsettings
    !     Calls:         -
    !     Arguments:     -
    !     Returns:       -
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff 2014
    !
    !-------------------------------------------------------------------c

    use m_filenames
    use m_optimization
    use m_atomproperties

    implicit none

    logical typelog,foundvariable,equals
    integer i,typeint,IOstatus,nwords,VarNamePlusEqLength,StringLength,wordLength
    real(8) typedble

    character*80 string
    character*80 VarName,word,VarValue,VarNamePlusEq
    character*80, allocatable:: words(:)
    character*1 tmpstring

    open(unit=1,file=trim(settingsfile))

    VarName="PARASTOOPT"
    foundvariable=.false.
    VarNamePlusEq=trim(VarName)//"="
    do
       read(1,'(a80)',IOSTAT=IOstatus) string
       !Currently not working:
       if (IOstatus.lt.0) then
          exit
       endif
       !if (string(1:9).eq."ENDOFFILE") then
       !   exit
       !endif

       call optionalarg(string,nwords)
       if (nwords.gt.0) then
        allocate(words(1:nwords))
        read(string,*) words

        if (words(1).ne.'#') then

          if (nwords.ge.1) then
             !First check if first word has an equals sign in it (eg,
             !VERBOSE=TRUE)
             word=words(1)
             wordLength=len_trim(word)
             equals=.false.
             do i=1,wordLength
               if (word(i:i).eq."=") then
                  equals=.true.
                  !print *,'equals in first word.'
                  exit
               endif
             enddo
             if (equals) then
               if (trim(word(1:i-1)).eq.trim(VarName)) then
                 !print *,trim(word(1:i-1)),'=',trim(VarName)
                 !print *,'i=',i,' wordLength=',wordLength
                 if (i.lt.wordLength) then
                   !print *,'..found it within first word..'
                   foundvariable=.true.
                   exit
                 else
                   if (nwords.ge.2) then
                     !print *,'..found it in second word..'
                     foundvariable=.true.
                     exit
                   endif
                 endif
               endif
             else
               if (trim(words(1)).eq.trim(VarName)) then
                 !print *,trim(words(1)),'=',trim(VarName)
                 if (nwords.ge.2) then
                   if (trim(words(2)).eq."=") then
                     if (nwords.ge.3) then
                       !print *,'..found it in third word..'
                       foundvariable=.true.
                       exit
                     endif
                   else
                     word=words(2)
                     wordLength=len_trim(word)
                     do i=1,wordLength
                       if (word(i:i).eq."=") then
                         equals=.true.
                         !print *,'equals in second word.'
                         exit
                       endif
                     enddo
                     if (equals) then
                       !print *,'..found it in second word..'
                       foundvariable=.true.
                       exit
                     endif
                   endif
                 endif
               endif
             endif
          endif

        endif
        deallocate(words)
       endif
    enddo
    read(1,*) freep(1:m3)  !cmin(maxspecies,maxspecies,maxspecies)
    read(1,*) freep(1+m3:2*m3) !cmax(maxspecies,maxspecies,maxspecies)
    read(1,*) freep(2*m3+1:2*m3+lm1*m1) !meamtau(0:lmax,maxspecies)
    read(1,*) freep(2*m3+lm1*m1+1:2*m3+lm1*m1+12*lm1*m2)
    !meamrhodecay(6,0:lmax,maxspecies,maxspecies)
    !meamrhodecay(6,0:lmax,maxspecies,maxspecies)
    read(1,*) freep(2*m3+lm1*m1+12*lm1*m2+1:2*m3+m1 &
        +lm1*m1+12*lm1*m2)   !meame0(maxspecies)
    read(1,*) freep(2*m3+m1+lm1*m1+12*lm1*m2+1:2*m3+ &
        2*m1+lm1*m1+12*lm1*m2) !meamrho0(maxspecies)
    read(1,*) freep(2*m3+2*m1+lm1*m1+12*lm1*m2+1:2*m3+ &
        3*m1+lm1*m1+12*lm1*m2) !meamemb3
    read(1,*) freep(2*m3+3*m1+lm1*m1+12*lm1*m2+1:2*m3+ &
        4*m1+lm1*m1+12*lm1*m2) !meamemb4
    read(1,*) freep(2*m3+(4+lm1)*m1+12*lm1*m2+1:m4)
    !pairpotparameter(16,maxspecies,maxspecies)
    read(1,*) freep(m4+1:m4+m2) !rs(maxspecies,maxspecies)
    read(1,*) freep(m4+m2+1:m4+2*m2) !rc(maxspecies,maxspecies)
    read(1,*) freep(m4+2*m2+1:m4+2*m2+m1) !enconst(maxspecies)
    close(1)

end subroutine readinoptparasonebyone
subroutine readmeamparam(filename)

    !---------------------------------------------------------------c
    !
    !     Read in a set of MEAM parameters from a file into the
    !     variables: cmin, etc. The name of the file is stored in
    !     the variable, filename.
    !
    !     Notes: lmaxPot is read in rather than lmax because lmax
    !     will already have been read in from the potential
    !     file (by extractlmax called from MEAMfit), but might have
    !     been changed by the settings file. E.g, we might want to
    !     optimize a MEAM potential but the input potential is EAM,
    !     in which case we want to keep lmax=3, but use a different
    !     varaible, lmaxPot, to denote which angular momenta of
    !     data we should read in from the potential file.
    !
    !     Called by:     program MEAMfit,initializemeam,meamenergy
    !     Calls:         -
    !     Arguments:     filename
    !     Returns:       cmin,cmax,meamtau,
    !                 meamrhodecay,meame0,meamrho0,
    !                 meamemb3,meamemb4,pairpotparameter,
    !                 rs,rc
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Ian Duff 2014
    !
    !---------------------------------------------------------------c

    use m_filenames
    use m_meamparameters
    use m_optimization
    use m_atomproperties

    implicit none

    integer i,j,l,lmaxPot
    real(8), allocatable:: meamtau_dummy(:,:)

    character*80 filename

    open(unit=1,file=trim(filename),status='old')
    read(1,*) lmaxPot
    !print *,'lmaxPot=',lmaxPot
    read(1,*) cmin !species 3x
    !print *,'cmin=',cmin
    read(1,*) cmax !species 3x
    !print *,'cmax=',cmax
    read(1,*) envdepmeamt !environ dep meamt's?
    !print *,'envdepmeamt=',envdepmeamt
    !meamtau is stored in the potential files with indices reversed (I.e., with
    !species first), for legacy reasons. Therefore first read into dummy array
    !before transfering to meamtau.
      !print *,'maxspecies=',maxspecies,' lmaxPot=',lmaxPot
    allocate(meamtau_dummy(1:maxspecies,0:lmaxPot))
    read(1,*) meamtau_dummy(1:maxspecies,0:lmaxPot)
    do i=1,maxspecies
       do l=0,lmaxPot
          meamtau(l,i)=meamtau_dummy(i,l)
       enddo
    enddo
    deallocate(meamtau_dummy)
    !print *,'meamtau=',meamtau
    !stop
    read(1,*) thiaccptindepndt
    !print *,'thiaccptindepnt=',thiaccptindepndt
    if (thiaccptindepndt.eqv..true.) then
        do i=1,maxspecies
            read(1,*) typethi(1,i)
            read(1,*) meamrhodecay(1:12,0:lmaxPot,1,i)
            !print *,'meamrhodecay=',meamrhodecay
        enddo
    else
        do i=1,maxspecies
            do j=1,maxspecies
                read(1,*) typethi(i,j)
                read(1,*) meamrhodecay(1:12,0:lmaxPot,i,j)
            enddo
        enddo
    endif
    read(1,*) embfunctype
    read(1,*) meame0 !species
    !print *,'meame0=',meame0
    read(1,*) meamrho0 !species
    !print *,'meamrho0=',meamrho0
    read(1,*) meamemb3 !species
    read(1,*) meamemb4 !species
    do i=1,maxspecies
        do j=1,maxspecies
            if (j.ge.i) then
                read(1,*) typepairpot(i,j)
                read(1,*) pairpotparameter(1:32,i,j)
    !print *,'pairpotparameter=',pairpotparameter
            endif
        enddo
    enddo
    read(1,*) rs !species 2x
    !print *,'rs=',rs
    read(1,*) rc !species 2
    !print *,'rc=',rc
    read(1,*) enconst
    !print *,'enconst=',enconst

    close(unit=1)

end subroutine readmeamparam
subroutine readoutcarenergy(filein,truedata_readin,freeEn,nconf,vasprun,currentposn)

    !----------------------------------------------------------------------c
    !
    !     Scan through the outcar file, 'filein', and locate the 'nconf'th
    !     true energy, which is returned as, 'truedata_readin'.
    !
    !     Called by:     readdatapoints
    !     Calls:         optionalarg
    !     Arguments:     filein
    !     Returns:       truedata_readin
    !     Files read:    filein
    !     Files written: -
    !
    !     Andy Duff, Jan 2008
    !
    !----------------------------------------------------------------------c
    implicit none

    logical freeEn,vasprun,skipnextposn
    integer i,j,ispc,it1,marker,nargs,nconf,iconf,wordLength,lineNmr,currentposn
    character*(*) filein
    character*9 newstrng,newstrng2
    character*80 line,strng1,strng2,aux(6),lineprev(30)
    character*80, allocatable:: words(:)
    real(8) truedata_readin

    ! print *,'Entering readoutcarenergy subroutine'

    if (currentposn.eq.0) then
       close(1)
       open(unit=1,file=filein)
       rewind(unit=1)
       iconf=1
    else
       iconf=currentposn
    endif

    do
       if (vasprun.eqv..false.) then 
          !Read in a line from the file and determine how many words
          !separated
          !by spaces there are in the file
          read(1,'(a80)') line
          call optionalarg(line,nargs)
          !Scan through the line to try and find the true energy
          allocate(words(nargs))
          read(line,*) words(1:nargs)
          if (freeEn.eqv..false.) then
             if (line(1:26).eq.'  energy  without entropy=') then
                 if (iconf.eq.nconf) then
                    read(line,*) aux(1:6),truedata_readin
                    exit
                 else
                    iconf=iconf+1
                 endif
             endif
          else
             if (line(1:25).eq.'  free  energy   TOTEN  =') then
                 if (iconf.eq.nconf) then
                    read(line,*) aux(1:4),truedata_readin
                    exit
                 else
                    iconf=iconf+1
                 endif
             endif
          endif
          deallocate(words)
       else

          skipnextposn=.false.
          do
            read(1,'(a80)') line
            if (line(1:10).eq.'<modeling>') then
               !This occurs at the start of a new vasprun.xml; which contains
               !the same initial set of positions twice
               skipnextposn=.true.
               ! print *,'Found <modeling>'
            endif
            newstrng=line(17:25)
            newstrng2=line(18:26)
            if ( (newstrng(1:9).eq.'positions') .or. &
                 (newstrng2(1:9).eq.'positions') ) then
              ! print *,'found positions'
              if (skipnextposn.eqv..false.) then
                 ! print *,'iconf=',iconf,', nconf=',nconf
                 if (iconf.eq.nconf) then
                    exit
                 endif
                 iconf=iconf+1
              else
                 skipnextposn=.false.
              endif
            endif
            call storelines(lineprev)
            lineprev(1)=line
          enddo

          !do i=1,30
          !print *,'lineprev(',i,')=',lineprev(i)
          !enddo

          !     <i name="e_0_energy">   -635.49695233</i>                                  
          if (freeEn.eqv..false.) then
             !Search for the word: e_0_energy in columns 14-23 in the past 30
             !lines (this expression can change positions based on VASP
             !settings), so can't just trace back a fixed no. of lines.
             do i=1,30
               line=lineprev(i)
               if (line(14:23).eq."e_0_energy") then
                 lineNmr=i
                 exit
               endif
             enddo
             line=lineprev(lineNmr) !17 before
             read(line,*) aux(1:2),strng2
             call keepNumber(strng2)
             strng2=trim(strng2)
             !wordLength=len_trim(strng2)
             !strng2=strng2(1:wordLength-1)
             read(strng2,*) truedata_readin
             ! print *,'energy=',truedata_readin
             ! stop
             exit
          else
             !Search for the word: e_fr_energy in columns 14-24 in the past 30
             !lines (this expression can change positions based on VASP
             !settings), so can't just trace back a fixed no. of lines.
             do i=1,30
               line=lineprev(i)
               if (line(14:24).eq."e_fr_energy") then
                 lineNmr=i
                 exit
               endif
             enddo 
             line=lineprev(lineNmr) !19 before
             read(line,*) aux(1:2),strng2
             call keepNumber(strng2)
             strng2=trim(strng2)
             !wordLength=len_trim(strng2) 
             !strng2=strng2(1:wordLength-1)
             read(strng2,*) truedata_readin
             ! print *,'read outcar energy, =',truedata_readin
             ! stop
             exit
          endif

       endif

    enddo

end subroutine readoutcarenergy


      subroutine storelines(lineprev)

      integer i
      character*80 lineprev(30)

      !Store 18 previously read lines (so we can back-trace to read in
      !the energy upon finding the atomic positions)

      do i=30,2,-1
         lineprev(i)=lineprev(i-1)
      enddo

      end subroutine storelines

subroutine readoutcarforces(filein,nforces,truedata_readin,istr,nconf,vasprun,currentposn)

    !----------------------------------------------------------------------c
    !
    !     Scan through the outcar file, 'filein', and locate the 'nconf'th
    !     set of true
    !     forces. Store these in the array, 'truedata_readin' (size
    !     3*nforces). If a row contains four columns, then the corresponding
    !     force is to be fit to ( optforce(i,istr)=.true. )
    !
    !     Called by:     readdatapoints
    !     Calls:         optionalarg
    !     Arguments:     filein,nforces
    !     Returns:       truedata_readin
    !     Files read:    filein
    !     Files written: -
    !
    !     Andy Duff, Jan 2008
    !
    !----------------------------------------------------------------------c

    use m_datapoints

    implicit none

    logical noForcesFlagged,vasprun
    integer i,nforces,nargs,istr,nconf,iconf,wordLength,currentposn
    character*1 snglchar
    character*80 strng1,strng2
    character*80 filein,line,lineprev(20)
    character*9 newstrng,newstrng2
    real(8) aux(1:4)
    real(8) truedata_readin(*)

    noForcesFlagged=.true.
    !print *,'startprevposn=',startprevposn
    if (currentposn.eq.0) then
       close(1)
       open(unit=1,file=filein)
       rewind(unit=1)
       iconf=1
    else
       iconf=currentposn
    endif

    do

       if (vasprun.eqv..false.) then

          !Read in a line from the file and determine how many words
          !separated by spaces there are in the file
          read(1,'(a80)') line
          call optionalarg(line,nargs)
          !         print *,'line=',line,', nargs=',nargs
          !Scan through the line to try and find the true forces
          !         allocate(words(nargs))
          !         read(line,*) words(1:nargs)
          !         print *,'words(1:',nargs,')=',words(1:nargs)

          !Insert code here to find the 'nconf'th set of forces

          if (line(1:9).eq.' POSITION') then
              if (iconf.lt.nconf) then
                 iconf=iconf+1
              else
                 !If 'TOTAL-FORCE' has been found, then read in the next
                 !lines as the forces
                 read(1,*)
                 !            stop
                 do i=1,nforces
                     !               print *,i,'/',nforces
                     !First determine if the line has a T at the end (to
                     !denote that this is a force to fit to)
                     read(1,*) line
                     if ((line(1:1).eq.'T').or.(line(1:1).eq.'t')) then
                         optforce(i,istr)=.true.
                         noForcesFlagged=.false.
                         backspace(unit=1)
                         read(1,*) snglchar,aux(1:3), &
                             truedata_readin(3*i-2:3*i)
                     else
                         optforce(i,istr)=.false.
                         backspace(unit=1)
                         read(1,*) aux(1:3), &
                             truedata_readin(3*i-2:3*i)
                     endif
                     !               call optionalarg(line,nargs)
                     !               print *,aux(1:3),truedata_readin(3*i-2:3*i-1)
                     !               print *,line,nargs
                     !               read(line,*) truedata_readin(3*i)
                     !               if (nargs.eq.2) then
                     !                  optforce(i,istr)=.true.
                     !                  print *,'fitting to force:',
                     !     +                truedata_readin(3*i-2:3*i)
                     !                  stop
                     !               else
                     !                  optforce(i,istr)=.false.
                     !               endif
                     !               stop
                     !               print *,truedata_readin(3*i-2:3*i)
                     !               stop
                     !            print *,'found force ',i,' (of ',nforces,') in OUTCAR file'
                 enddo
                 exit
              endif
          endif
          !         deallocate(words)
       else

          do
            read(1,'(a80)') line
            newstrng=line(17:25)
            newstrng2=line(18:26)
            if ( (newstrng(1:6).eq.'forces') .or. &
                 (newstrng2(1:6).eq.'forces') ) then
               if (iconf.eq.nconf) then
                  exit
               endif
               iconf=iconf+1
            endif
            call storelines(lineprev)
            lineprev(1)=line
          enddo
          !print *,iconf,'/',nconf
          do i=1,nforces
             !of form "   <v>       0.00000213      0.00000295      -0.00000143</v>"
             !print *,i,'/',nforces 
             read(1,*) line
             if ((line(1:1).eq.'T').or.(line(1:1).eq.'t')) then
                 optforce(i,istr)=.true.
                 noForcesFlagged=.false.
                 backspace(unit=1)
                 read(1,*) snglchar,strng1,truedata_readin(3*i-2:3*i-1),strng2
             else
                 optforce(i,istr)=.false.
                 backspace(unit=1)
                 read(1,*) strng1,truedata_readin(3*i-2:3*i-1),strng2
             endif
             call keepNumber(strng2)
             strng2=trim(strng2)
             !wordLength=len_trim(strng2)
             !strng2=strng2(1:wordLength-1)
             read(strng2,*) truedata_readin(3*i)
             !print *,'force(',i,')=',truedata_readin(3*i-2:3*i)

          enddo
          exit 

       endif

    enddo

    if (noForcesFlagged) then
        write(*,'(A29,I2)') ' Fitting all forces for file ',istr
        do i=1,nforces
           optforce(i,istr)=.true.
        enddo
    else
        write(*,'(A47,I2)') 'Fitting only manually selected forces for file ',istr
    endif

end subroutine readoutcarforces
subroutine readposcar(filein,nconf,vasprun,startprevposn)

    !-----------------------------------------------------------c
    !
    !     Reads VASP type vasprun.xml files or alternatively
    !     POSCAR or CONTCAR files and finds the
    !     'nconfig'th ionic configuration in the file.
    !     Converts atomic coordinates to direct coordinates if
    !     necessary and then makes sure all atoms are in the
    !     central unit cell. Finally atomic coordinates are
    !     converted to cartesian coordinates
    !
    !     Called by:     initializestruc
    !     Calls:         optionalarg,matinv,matmul
    !     Arguments:     filein
    !     Returns:       tr,nspecies,nat,natoms,z,zz,coordinates,
    !                    coords_old
    !     Files read:    filein
    !
    !     Marcel Sluiter, Feb 8 2006
    !
    !-----------------------------------------------------------c

    use m_atomproperties
    use m_poscar

    implicit none

    logical num,startprevposn,vasprun,skipnextposn
    logical, parameter:: force_atoms_into_cell=.false. !If true, the
    !atoms in the poscar file are translated
    !by the basis vectors so that they lie
    !within the unit cell defined by these
    !basis vectors
    integer iconf,i,j,it1,nconf,wordLength
    real(8) aux1(3),aux2(3),scale
    character*(*) filein
    character*9 newstrng,newstrng2
    character*80 strng1,strng2
    character*80 line,zspecies
    character*80, allocatable:: words(:)

    if (startprevposn.eqv..false.) then
       close(1)
       open(unit=1,file=filein)
       !print *,'reading from scratch from file:',trim(filein)
       if (nconf.gt.1) then
          if (vasprun.eqv..true.) then
             !Skip past the steps we are not interested in.
             iconf=1
             skipnextposn=.false.
             do
               read(1,'(a80)') line
               if (line(1:10).eq.'<modeling>') then
                  !print *,'found string <modeling>'
                  !This occurs at the start of a new vasprun.xml; which contains
                  !the same initial set of positions twice
                  skipnextposn=.true.
                  call readAtomPropsFromVasprun
                  !print *,'nspecies=',nspecies
                  !print *,'natoms=',natoms
                  !print *,'z=',z
                  !print *,'nat=',nat
                  !print *,'zz=',zz
                  !print *,'tr=',tr
               endif
               newstrng=line(17:25)
               newstrng2=line(18:26)
               if ( (newstrng(1:9).eq.'positions') .or. &
                    (newstrng2(1:9).eq.'positions') ) then
                 if (skipnextposn.eqv..false.) then
                    if (iconf.eq.(nconf-1)) then
                       exit
                    endif
                    iconf=iconf+1
                 else
                    skipnextposn=.false.
                 endif
               endif
             enddo
          else
             !Skip past the ionic steps we are not interested in
             do iconf=1,nconf-1
                read(1,'(a80)') zspecies
                read(1,*)
                read(1,*)
                read(1,*)
                read(1,*)
                read(1,'(a80)') line
                call optionalarg(line,nspecies)
                call checkifnumber(line,num) !Check if line has species types
                if (num.eqv..false.) then
                   read(1,'(a80)') line
                endif
                read(line,*) nat(1:nspecies)
                natoms = sum( nat(1:nspecies) )
                read(1,*) line            !cartesian or direct
                do i=1,natoms
                    read(1,*)
                enddo
             enddo
          endif
       endif
    endif

    if (vasprun.eqv..false.) then
       read(1,'(a80)') zspecies   !This is a VASP comment line BUT we are
       !going to use it to indicate the Z of each atomic species.
       read(1,*) scale
       read(1,*) tr(1:3,1)       !a lattice parameter
       read(1,*) tr(1:3,2)       !b lattice parameter
       read(1,*) tr(1:3,3)       !c lattice parameter
       tr=tr*scale
       read(1,'(a80)') line      !Numbers of atoms of each species (or species types)
       call optionalarg(line,nspecies)
       if(allocated(nat)) deallocate(nat)
       if(allocated(z)) deallocate(z)
       if(allocated(coordinates)) deallocate(coordinates)
       if(allocated(zz)) deallocate(zz)
       allocate(nat(nspecies),z(nspecies))
       call checkifnumber(line,num) !Check if line has species types
       if (num.eqv..false.) then
           allocate(words(nspecies))
           read(line,*) words(1:nspecies)
           call elementToZ(words)
           deallocate(words)
           !print *,'species as read from line above atom numbers:',z(1:nspecies)
           read(1,'(a80)') line      !Numbers of atoms of each species
       else
           call checkifnumber(zspecies,num)
           if (num.eqv..false.) then
              allocate(words(nspecies))
              read(zspecies,*) words(1:nspecies)
              call elementToZ(words)
              deallocate(words)
           else
              read(zspecies,*) z(1:nspecies) !z=14
           endif
       endif
       read(line,*) nat(1:nspecies)
       !z(1:nspecies)
       natoms = sum( nat(1:nspecies) )
    
       !Fill zz array:
       allocate(zz(natoms))
       it1=0
       do i=1,nspecies
           do j=1,nat(i)
               it1=it1+1
               zz(it1)=z(i)
               !print *,'zz(',it1,')=',zz(it1),'=z(',i,')=',z(i)
           enddo
       enddo
       stop
       read(1,*) line            !cartesian or direct
       allocate(coordinates(3,natoms))
       do i=1,natoms
           read(1,*) coordinates(1:3,i)
       enddo
    
       !To ensure all atoms are in the central unit cell, we first go to
       !direct coordinates. If cartesian coordinates, convert:
       if(line(1:1).eq.'C'.or.line(1:1).eq.'c') then
           coordinates=coordinates*scale
           call matinv(tr,3,3,trinv)
           do i=1,natoms
               aux1 = coordinates(1:3,i)
               aux2 = matmul(trinv,aux1)
               coordinates(1:3,i) = aux2
           enddo
       endif
    
    else

       !insert code to read in from vasprun the following:
       !nspecies, natoms, z, nat, zz, tr, trinv and
       !coordinates (in direct coords)

          !Skip past the steps we are not interested in.
          skipnextposn=.false.
          !print *,'looking for <modeling>'
          do
            read(1,'(a80)') line
            if (line(1:10).eq.'<modeling>') then
               !print *,'found <modeling>'
               !This occurs at the start of a new vasprun.xml; which contains
               !the same initial set of positions twice
               skipnextposn=.true.
               !2015 July: Andy added---
               if(allocated(nat)) deallocate(nat)
               if(allocated(z)) deallocate(z)
               if(allocated(zz)) deallocate(zz)
               !---
               call readAtomPropsFromVasprun
               !print *,'returned!'
               !print *,'nspecies=',nspecies
               !print *,'natoms=',natoms
               !print *,'z=',z
               !print *,'nat=',nat
               !print *,'zz=',zz
               !print *,'tr=',tr
            endif
            newstrng=line(17:25)
            newstrng2=line(18:26)
            if ( (newstrng(1:9).eq.'positions') .or. &
                 (newstrng2(1:9).eq.'positions') ) then
              !print *,'found positions'
              if (skipnextposn.eqv..false.) then
                 exit
              else
                 !print *,'skipping next posn...'
                 skipnextposn=.false.
              endif
            endif
          enddo
          !print *,'left first loop'
          if(allocated(coordinates)) deallocate(coordinates)
          allocate(coordinates(3,natoms))
          !print *,'allocated arrays, preping to read coords'
          do i=1,natoms
             read(1,*) strng1,coordinates(1,i),coordinates(2,i),strng2
             call keepNumber(strng2)
             strng2=trim(strng2)
             !wordLength=len_trim(strng2)
             !strng2=strng2(1:wordLength-1)
             read(strng2,*) coordinates(3,i)
             !print *,'i=',i,' coords=',coordinates(1:3,i)
          enddo
          !stop
          !print *,'read coords'

    endif

    if (force_atoms_into_cell) then
        !Use modulo functions to get all coordinates in 0-1 range
        do i=1,natoms
            do j=1,3
                coordinates(j,i)=modulo(coordinates(j,i),1d0) !modulo
                !gives sign of 2nd argument
            enddo
        enddo
    endif
    !Convert direct coordinates to cartesian
    do i=1,natoms
        aux1 = coordinates(1:3,i)
        aux2 = matmul(tr,aux1)
        coordinates(1:3,i) = aux2
    enddo
    !Store the coordinates in 'coords_old' for the case where we are
    !calculating forces. We need to keep the old coordinates and no.
    !of atoms in unit cell handy for when we perform the second scan
    !across trial unit cells to establish which atoms satisfy
    !p_rmax1 < dist < p_rmax2
    if(allocated(coords_old)) deallocate(coords_old)
    allocate(coords_old(3,natoms))
    natoms_old=natoms
    do i=1,natoms
        coords_old(1:3,i)=coordinates(1:3,i)
    enddo

end subroutine readposcar

subroutine checkifnumber(line,num)

    implicit none

    logical num
    integer i,length
    character*(*) line

    length=len(line)
    !print *,'line=',line,', length=',length
    if (length.eq.0) then
       print *,'ERROR: length=0 in checkifnumber, STOPPING.'
       stop
    endif
   
    num=.false.
    i=1
    do
       if (line(i:i).eq.' ') then
          if (i.eq.length) then
             exit
          endif
          i=i+1
       else
          if (line(i:i).eq.'1') then
             num=.true.
          elseif (line(i:i).eq.'2') then
             num=.true.
          elseif (line(i:i).eq.'3') then
             num=.true.
          elseif (line(i:i).eq.'4') then
             num=.true.
          elseif (line(i:i).eq.'5') then
             num=.true.
          elseif (line(i:i).eq.'6') then
             num=.true.
          elseif (line(i:i).eq.'7') then
             num=.true.
          elseif (line(i:i).eq.'8') then
             num=.true.
          elseif (line(i:i).eq.'9') then
             num=.true.
          elseif (line(i:i).eq.'0') then
             num=.true.
          endif
          exit
       endif
    enddo

endsubroutine checkifnumber


subroutine elementToZ(words)

    use m_atomproperties
    use m_poscar

    implicit none

    integer i,j
    character*80 words(nspecies)

    do i=1,nspecies
        do j=1,88
           if (words(i).eq.element(j)) then
              z(i)=j
           endif
        enddo
    enddo

end subroutine elementToZ

subroutine readsettings(nfreep,iseed)

    !-------------------------------------------------------------------c
    !
    !     Reads in the settings from the settings file.
    !
    !     Called by:     program MEAMfit
    !     Calls:         -
    !     Arguments:     settingsfile,maxspecies,lmax
    !     Returns:       nsteps,freep
    !     Files read:    settingsfile
    !     Files written: -
    !
    !     Andrew Duff 2006
    !
    !-------------------------------------------------------------------c

    use m_filenames
    use m_datapoints
    use m_atomproperties
    use m_generalinfo
    use m_meamparameters
    use m_geometry
    use m_optimization
    use m_plotfiles

    implicit none

    logical onebyone,exist,foundvariable,opt_enconst,opt_meamtau
    integer i,j,k,nfreep,iseed,nPara,ncutoffs,&
            ipara,paraCnt,l,nembterms,ncutoffDens_overall,nterms_pp_overall
    integer, allocatable:: ncutoffDens(:,:),ncutoffPairpot(:,:), &
            nterms_pp(:,:)
    character*1 string2,string3
    character*5 writecheck
    character*80 searchfor,found
    character*80 string

    !Check there is a settings file. If not, write one and stop.
    open(unit=1,file=trim(settingsfile))
    if (settingsfileexist.eqv..false.) then
       print *
       print *,'WARNING: settings file does not exist.'
       upplimoptfunc=10.d0
       maxoptfuncallowed=0.1d0
       optimizestruc=.false.
       nsteps=100000
       lookuptables=.false.
       onebyone=.true.
       write(*,*) 'Writing template settings file.'
       write(1,*) 'TYPE=EAM'
       write(1,*) '# CUTOFF_MAX='
       write(1,*) 'NTERMS=2'
       write(1,*) 'NTERMS_EMB=3'
       stop
    else

       !Search through settings file tags (VERBOSE, etc) and place values in the
       !corresponding variables
       allocate( typethi(maxspecies,maxspecies), &
                 typepairpot(maxspecies,maxspecies) )
       readParasOnebyone=.false.
       printoptparas=.false.

       searchfor="VERBOSE" !Extra information written to screen
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) verbose
          write(*,'(A9,L1,A21)') ' verbose=',verbose,' (from settings file)'
       else
          verbose=.false.
       endif

       searchfor="FASTFORCE" !Store separations for displaced configurations
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) fastForce
          write(*,'(A11,L1,A21)') ' FASTFORCE=',fastForce,' (from settings file)'
       else
          fastForce=.true.
       endif

       searchfor="USEREF" !Use reference structure to calculate energy differences
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) useRef
          write(*,'(A8,L1,A21)') ' USEREF=',useRef,' (from settings file)'
       else
          if (verbose) then
             print *,'No USEREF in settings file, using USEREF=true'
          endif
          useRef=.true.
       endif

       searchfor="SEED" !A seed can be specified for debugging purposes
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) iseed
          write(*,'(A6,I10,A21)') ' seed=',iseed,' (from settings file)'
       else
          print *,'No SEED in settings file, using time to seed random numbers.'
          iseed=0
       endif

       searchfor="MUT_PROB" !Mutation probability in the GA
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) probRand
          write(*,'(A10,F15.10,A21)') ' MUT_PROB=',probRand,' (from settings file)'
       else
          print *,'No MUT_PROB in settings file, using default (0.3).'
          probRand=0.3d0
         !if (verbose) then
         !   print *,'   (Probability of a mutation, I.e. random parameterization, for a'
         !   print *,'    given potential parameter/group of potential parameters)'
         !endif
       endif

       searchfor="PROBPOT1" !Probability that 'genes' of parent 1 are passed to offspring
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) probPot1
          write(*,'(A10,F15.10,A21)') ' PROBPOT1=',probPot1,' (from settings file)'
       else
          if (verbose) then
             print *,'No PROBPOT1 in settings file, using default (0.5).'
          endif
          probPot1=0.5d0
       endif

       searchfor="NOPTFUNCSTORE" !No. parents to be stored for GA (also the number of 
                                 !potentials to be saved)
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) noptfuncstore
          write(*,'(A15,I10,A21)') ' NOPTFUNCSTORE=',noptfuncstore,' (from settings file)'
       else
          print *,'No NOPTFUNCSTORE in settings file, using default (10).'
          noptfuncstore=10
       endif

       searchfor="OPTFUNCCG" !Optfunc must be less than this to move from random 
                             !sampling to CG optimization phase
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) upplimoptfunc
          write(*,'(A11,F10.5,A21)') ' OPTFUNCCG=',upplimoptfunc,' (from settings file)'
       else
          print *,'No OPTFUNCCG in settings file, using default (OPTFUNCCG=10).'
          upplimoptfunc=10
       endif

       searchfor="GENALGO" !Whether to use genetic algorithm having filled up
                           !optimization function table
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) genAlgo
          write(*,'(A9,F10.5,A21)') ' GENALGO=',genAlgo,' (from settings file)'
       else
          if (verbose) then
             print *,'No GENALGO in settings file, using default (GENALGO=TRUE).'
          endif
          genAlgo=.true.
       endif

       searchfor="OPTFUNC_ERR" !This quantity is used in the CG routine, and denotes
                             !the error in calculating the optimization
                             !function. Needs reducing for force-fit due to
                             !finite-difference approach
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) optfunc_err
          write(*,'(A13,F10.5,A21)') ' OPTFUNC_ERR=',optfunc_err,' (from settings file)'
       else
          if (forcefit.eqv..false.) then
             print *,'No OPTFUNC_ERR in settings file, using default for energy fit (OPTFUNC_ERR=10^-14).'
             optfunc_err=1.d-14
          else
             print *,'No OPTFUNC_ERR in settings file, using default for fit including forces (OPTFUNC_ERR=10^-4).'
             optfunc_err=1.d-4
          endif
       endif

       searchfor="OPTFUNCCG_GA" !Optfunc must be less than this for a spliced/mutated 
                                !offspring to be CG optimized
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) upplimoptfunc_GA
          if (verbose) then
             write(*,'(A14,F10.5,A21)') ' OPTFUNCCG_GA=',upplimoptfunc_GA,' (from settings file)'
          endif
       else
          write(*,'(A60)') ' No OPTFUNCCG_GA in settings file, using default (OPTFUNCCG='
          write(*,'(A13,F15.10,A2)') ' 10*OPTFUNCCG=',10d0*upplimoptfunc,').'
          upplimoptfunc_GA=10d0*upplimoptfunc
       endif

       searchfor="OPTDIFF" !The precision in the relative optimization function at which
                           !CG minimization is terminated
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) optdiff
          write(*,'(A9,E10.5,A21)') ' OPTDIFF=',optdiff,' (from settings file)'
       else
          print *,'No OPTDIFF in settings file, using default (=10^-10)'
          optdiff=1d-10
       endif

       searchfor="OPTACC" !When the optfuncs of the best and worst stored potentials are within
                          !this amount of one another, optimization will stop.
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) optAcc
          write(*,'(A8,F10.5,A21)') ' OPTACC=',optAcc,' (from settings file)'
       else
          print *,'No OPTACC in settings file, using default (=0.0005)'
          optAcc=0.0005d0
       endif

       searchfor="STOPTIME" !Optimization will stop after this length of time (although a given
                            !CG optimization will be allowed to finish first)
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) stopTime
          write(*,'(A10,I5,A21)') ' STOPTIME=',stopTime,' (from settings file)'
       else
          print *,'No STOPTIME in settings file, using default (=168 hours=1 week)'
          stopTime=168
       endif

       searchfor="MAXFUNCEVALS" !Maximum number of function evaluations before CG minimzation
                                !is terminated)
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) maxfuncevals
          if (verbose) then
             write(*,'(A14,I7,A21)') ' MAXFUNCEVALS=',maxfuncevals,' (from settings file)'
          endif
       else
          print *,'No MAXFUNCEVALS in settings file, using default (=2000)'
          maxfuncevals=2000
       endif

       searchfor="OPTFUNCSTP" !If the optfunc is smaller than this, optimization will stop
                              !(Andy: no practical use - remove this)
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) maxoptfuncallowed
          write(*,'(A12,F10.5,A21)') ' OPTFUNCSTP=',maxoptfuncallowed,' (from settings file)'
       else
          if (verbose) then
             print *,'No OPTFUNCSTP in settings file, using default (=0d0)'
          endif
          maxoptfuncallowed=0d0
       endif

       searchfor="OPTIMIZESTRUC" !If true, relax crystal structures according to potential during
                                 !optimization (Andy: do not use yet, needs further development/testing)
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) optimizestruc
          write(*,'(A13,L1,A21)') ' OPTIMIZESTRUC=',optimizestruc,' (from settings file)'
       else
          if (verbose) then
             print *,'No OPTIMIZESTRUC in settings file, using default (=.false.)'
          endif
          optimizestruc=.false.
       endif

       !Check if a potential is to be read in from a file
       if (readpotfile.eqv..true.) then
          print *,'POTFILEIN=',trim(startparameterfile),' (from command-line)'
       else
          searchfor="POTFILEIN"
          call getsettingspara(searchfor,startparameterfile,foundvariable)
          if (foundvariable) then
             print *,'POTFILEIN=',trim(startparameterfile),' (from settings)'
             readpotfile=.true.
          else
             if (verbose) then
                print *,'No POTFILEIN in settings file, not reading in potential from file'
             endif
             readpotfile=.false.
          endif
       endif

       searchfor="POTFILEOUT" !Produce an example potential file as a guide for user to enter in
                              !their own starting point for an optimization (Andy: in development)
       call getsettingspara(searchfor,startparameterfile,foundvariable)
       if (foundvariable) then
          print *,'POTFILEOUT=',trim(startparameterfile),' (from settings file)'
          if (readpotfile) then
             print *,'ERROR: you cannot both read in (POTFILEIN) and read out (POTFILEOUT)'
             print *,'a potential file, stopping.'
             stop
          endif
          writepotfile=.true.
       else
          if (verbose) then
             print *,'No POTFILEOUT in settings file, using default (=.false.)'
          endif
          writepotfile=.false.
       endif

       searchfor="FIXPOTIN" !Check if the input potential is to be fixed or allowed to optimize
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          if (readpotfile.eqv..false.) then
             print *,'ERROR: FIXPOTIN specified, yet no input parameter file provided,'
             print *,'stopping.'
             stop
          endif
          read(found,*) fixPotIn
          write(*,'(A10,L1,A21)') ' FIXPOTIN=',fixPotIn,' (from settings file)'
       else
          if (readpotfile) then
             print *,'No FIXPOTIN in settings file: allow input potential to optimize'
          endif
          fixPotIn=.false.
       endif

       !Check if a continuation job is to be performed (uses the potparas_best files and 
       !bestoptfuncs file as starting point)
       if (contjob.eqv..true.) then
          write(*,'(A6,L1,A20)') ' CONT=',contjob,' (from command line)'
       else
          searchfor="CONT"
          call getsettingspara(searchfor,found,foundvariable)
          if (foundvariable) then
             read(found,*) contjob
             write(*,'(A6,L1,A21)') ' CONT=',contjob,' (from settings file)'
             if (readpotfile.and.contjob) then
                print *,'ERROR: you cannot have POTFILEIN=T and CONT=T simultaneously, stopping.'
                stop
             endif
          else
             if (verbose) then
                print *,'No CONT in settings file, using default (=.false.)'
             endif
             contjob=.false.
          endif
       endif

       !Check if an optimization is to be performed (NOPT=TRUE if not)
       if (noOpt.eqv..true.) then
          !Check first if NOOPT has been specified already from the command-line
          nsteps=1
          startparas=1
          write(*,'(A31)') ' NOOPT=TRUE (from command-line)'
       else
          !Default values if NOOPT is not found:
          nsteps=2 !infinite no. steps
          startparas=1
          searchfor="NOOPT" !aka nsteps, which is either 0 for no opt or >0
                               !otherwise
          call getsettingspara(searchfor,found,foundvariable)
          if (foundvariable) then
             read(found,*) noOpt
             if (noOpt.eqv..true.) then
                write(*,'(A32)') ' NOOPT=TRUE (from settings file)'
                nsteps=1
                startparas=1
             endif
          endif
       endif

       !Read in maximum cut-off radius - necessary if optimizing; writing an example potential
       !file; or producing a separation histogram. Cut-off radii cannot take values greater than
       !this.
       searchfor="CUTOFF_MAX"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) cutoffMax
          write(*,'(A12,F10.5,A21)') ' CUTOFF_MAX=',cutoffMax,' (from settings file)'
       else
          if ((nsteps.gt.1).or.(startparas.eq.2).or.(writepotfile)) then
             print *,'No CUTOFF_MAX in settings file, MUST be specified,'
             print *,'STOPPING.'
             stop
          endif
          if ((startparas.eq.1).and.(readpotfile.eqv..false.)) then
             print *,'No CUTOFF_MAX in settings file. Please use this tag to'
             print *,'supply a maximum seperation to be plotted out to in the'
             print *,'sepnHistogram file, STOPPING.'
             stop
          endif
       endif

       !Read in minimum cut-off radius if provided. Cut-off radii cannot take values smaller 
       !than this
       searchfor="CUTOFF_MIN"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) cutoffMin
          if (cutoffMin.lt.cutoffMinLim) then
             print *,'ERROR: cutoffMin < ',cutoffMinLim,'. Please choose larger value.'
             print *,'(cutoff radii would enter the Ziegler-Biersak region otherwise) Stopping.'
             stop
          endif
          write(*,'(A12,F10.5,A21)') ' CUTOFF_MIN=',cutoffMin,' (from settings file)'
       else
          print *,'No CUTOFF_MIN in settings file, using default (CUTOFF_MIN=',cutoffMinLim,')'
          cutoffMin=cutoffMinLim
       endif

       !Read in potentials from r vs function tables. Andy: not yet implemented
       searchfor="LOOKUPTABLES"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) lookuptables
          write(*,'(A12,L1,A21)') ' LOOKUPTABLES=',lookuptables,' (from settings file)'
       else
          if (verbose) then
             print *,'No LOOKUPTABLES in settings file, using default (=.false.)'
          endif
          lookuptables=.false.
       endif

    endif

    if (settingsfileexist.eqv..true.) then

       !Read in type of potential (EAM or MEAM)
       searchfor="TYPE"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) string
          if (string(1:3).eq.'EAM') then
             print *,'TYPE=EAM (from settings file)'
             if ((readpotfile.eqv..true.).and.(lmax.gt.0)) then
                print *,'ERROR: potential parameter input file describes a MEAM potential,'
                print *,'conflicting with this setting, STOPPING.'
                stop
             endif
             lmax=0
          elseif (string(1:4).eq.'MEAM') then
             print *,'TYPE=MEAM (from settings file)'
             lmax=3
          else
             print *,'ERROR: TYPE not recognized in settings file'
             stop
          endif
       else
          if (readpotfile.eqv..false.) then
             print *,'No TYPE in settings file, using default (=EAM)'
             lmax=0
          endif
       endif

       !Initialize freep array. This is the array which tells the code which
       !parameters are to be optimized (and will be filled in this subroutine -
       !but possibly amended again later in the code)
       np=2*maxspecies*maxspecies*maxspecies+ &
           maxspecies*maxspecies*34+ &
           maxspecies*(6+lmax)+ &
           12*(lmax+1)*(maxspecies*maxspecies)
       if (.not.allocated(freep)) then
           allocate(freep(np))
       endif
       freep=0
       call initParaLimits !Set up variables useful for navigating the freep() and p() arrays

       !Upper and lower bounds for random generation of potentials are usually
       !set in-code - the following option allows explicit specification of
       !these limits from the settings file
       searchfor="RANDPARASMANUAL"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) readparasfromsettings
          write(*,'(A17,L1,A21)') ' RANDPARASMANUAL=',readparasfromsettings,' (from settings file)'
       else
          if (verbose) then
             print *,'No RANDPARASMANUAL in settings file, using default (=.FALSE.)'
          endif
          readparasfromsettings=.false.
       endif

       !Option for using environment dependent meam_t's
       if (lmax.gt.0) then
          searchfor="ENVDEPMEAMT"
          call getsettingspara(searchfor,found,foundvariable)
          if (foundvariable) then
             read(found,*) envdepmeamt
             write(*,'(A13,L1,A21)') ' ENVDEPMEAMT=',envdepmeamt,' (from settings file)'
          else
             if (verbose) then
                print *,'No ENVDEPMEAMT in settings file, using default (=.FALSE.)'
             endif
             envdepmeamt=.false.
          endif
       endif

       !Option for allowing electron densities to depend not just on the species
       !of the atom from which it is originating but also on the species of atom
       !which it is incident on (set THIACCPTINDEPNDT=.false. if this is required).
       !Andy: please do not use yet - should work but not extensively tested.
       searchfor="THIACCPTINDEPNDT"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) thiaccptindepndt
          write(*,'(A18,L1,A21)') ' THIACCPTINDEPNDT=',thiaccptindepndt,' (from settings file)'
       else
          if (verbose) then
             print *,'No THIACCPTINDEPNDT in settings file, using default (=.TRUE.)'
          endif
          thiaccptindepndt=.true.
       endif

       !Analytic form of electron density (Andy: curently only one available but
       !future developments will enable multiple types - also for TYPE_EMB, TYPE_PP)
       searchfor="TYPE_DENS"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) typethi
          write(*,'(A11,I1,A21)') ' TYPE_DENS=',typethi,' (from settings file)'
       else
          if (verbose) then
             print *,'No TYPE_DENS in settings file, using default (=1)'
          endif
          typethi=1
       endif

       !Analytic form for embedding function
       searchfor="TYPE_EMB"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) embfunctype
          write(*,'(A10,I1,A21)') ' TYPE_EMB=',embfunctype,' (from settings file)'
       else
          if (verbose) then
             print *,'No TYPE_EMB in settings file, using default (=1)'
          endif
          embfunctype=1
       endif

       !Analytic form for pair-potential
       searchfor="TYPE_PP"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) typepairpot
          write(*,'(A9,I1,A21)') ' TYPE_PP=',typepairpot,' (from settings file)'
       else
          if (verbose) then
             print *,'No TYPE_PP in settings file, using default (=2)'
          endif
          typepairpot=2
       endif

       !Maximum electron density to be used in the LAMMPS and Camelion output
       !files. Can either be specified or will be determined based on the
       !maximum value encountered during the fitting.
       searchfor="RHOMAX_LAMMPS"
       call getsettingspara(searchfor,found,foundvariable)
       if (foundvariable) then
          read(found,*) rhomax
          write(*,'(A15,F20.10,A21)') ' RHOMAX_LAMMPS=',rhomax,' (from settings file)'
          autorhomax=.false.
       else
          if (verbose) then
             print *,'No RHOMAX_LAMMPS in settings file, will create'
             print *,'automatically based on maximum rho encountered from'
             print *,'fitting database.'
          endif
          autorhomax=.true.
       endif

      !searchfor="SPLN_NVALS"
      !call getsettingspara(searchfor,found,foundvariable)
      !if (foundvariable) then
      !   read(found,*) splnNvals
      !   write(*,'(A12,I10,A21)') ' SPLN_NVALS=',splnNvals,' (from settings file)'
      !else
      !   if (verbose) then
      !      print *,'No SPLN_NVALS in settings file, using default (=100000),'
      !   endif
      !   splnNvals=100000
      !endif

       if ((nsteps.gt.1).or.(startparas.eq.2)) then

        !See if the parameters to optimized are to be given explicitly to the
        !code (PARASTOOPT={...}, where ... are a list of 0's 1's and 2's in the
        !same order as the parameters are given in the potparas_best files. An
        !example of such a list can be requested by entering PARASTOOPT=write).
        !Andy: I do not expect people will use this, as the code has been
        !designed to avoid it - but it might be useful under specific
        !circumstances)
        searchfor="PARASTOOPT"
        call getsettingspara(searchfor,found,foundvariable)
        if (foundvariable) then
           if (fixPotIn.eqv..true.) then
              print *,'ERROR: Use of FIXPOTIN inconsistent with use of PARASTOOPT.'
              print *,'Please remove one or the other and restart. Stopping.'
              stop
           endif
           writecheck=found(1:5)
           if ((writecheck.eq."write").or.(writecheck.eq."WRITE").or.(writecheck.eq."Write")) then
              print *,'Found PARASTOOPT: Write flag detected. Outputing example input for this tag.'
              printoptparas=.true.
           else
              print *,'Found PARASTOOPT (reading optimization parameters in one-by-one)'
              call readinoptparasonebyone
           endif
           readParasOnebyone=.true.
           startparas=2
        else
           if (verbose) then
              print *,'No PARASTOOPT, searching for NTERM tags'
           endif
           readParasOnebyone=.false.
        endif
        allocate(ncutoffDens(maxspecies,maxspecies), &
                 ncutoffPairpot(maxspecies,maxspecies), &
                 nterms_pp(maxspecies,maxspecies))

        !If PARASTOOPT is not specified, instead look for NTERMS flags. NTERMS
        !specifies the number of terms to be used in the pairwise functions
        ncutoffs=0
        searchfor="NTERMS"
        call getsettingspara(searchfor,found,foundvariable)
        if (foundvariable) then
           if (readParasOnebyone) then
             print *,'ERROR: NTERMS tag conflicts with PARASTOOPT. Please remove'
             print *,'one and restart, stopping.'
             stop
           endif
           read(found,*) ncutoffs
           if (verbose) then
              write(*,'(A8,I2,A21)') ' NTERMS=',ncutoffs,' (from settings file)'
           endif
           startparas=2
        endif

        !NTERMS can instead be specified for the electron densities and
        !pair-potentials separately using NTERMS_DENS and NTERMS_PP
        ncutoffDens_overall=0
        searchfor="NTERMS_DENS"
        call getsettingspara(searchfor,found,foundvariable)
        if (foundvariable) then
           if (readParasOnebyone) then
             print *,'ERROR: NTERMS_DENS tag conflicts with PARASTOOPT. Please'
             print *,'remove one and restart, stopping.'
             stop
           endif
           read(found,*) ncutoffDens_overall
           if (verbose) then
              write(*,'(A13,I1,A21)') &
             ' NTERMS_DENS=',ncutoffDens_overall,' (from settings file)'
           endif
           startparas=2
        endif
        paraCnt=2*m3+lm1*m1+1
        do i=1,maxspecies
           do j=1,maxspecies
              if ((i.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                 write(string2,'(I1)') i
                 write(string3,'(I1)') j
                 !NTERMS_DENS tags can also be written as NTERMS_DENS(1,1), 
                 !or NTERMS_DENS(1,2) for example, which for a 2-species system
                 !would provide the numbers of terms in the electron densities
                 !for species 1 and 2 separately (the first index would be used
                 !in the case that THIACCPTINDEPNDT=.FALSE.)
                 searchfor="NTERMS_DENS("//string2//";"//string3//")"
                 call getsettingspara(searchfor,found,foundvariable)
                 if (foundvariable) then
                    if (readParasOnebyone) then
                      write(*,'(A8,A16,A31)') ' ERROR: ',trim(searchfor),' tag conflicts with PARASTOOPT.'
                      print *,'Please remove one and restart, stopping.'
                      stop
                    endif
                    read(found,*) ncutoffDens(i,j)
                    write(*,'(A13,I1,A1,I1,A2,I2,A21)') &
          ' NTERMS_DENS(',i,';',j,')=',ncutoffDens(i,j),' (from settings file)'
                    startparas=2
                 else
                    if ((readParasOnebyone.eqv..false.).and.(startparas.eq.2)) then
                       if (ncutoffDens_overall.gt.0) then
                          !If an NTERMS_DENS(i,j) has not been explicitly
                          !specified, use the 'master' NTERMS_DENS value.
                          write(*,'(A13,I1,A1,I1,A16,I2,A1)') &
                ' NTERMS_DENS(',i,';',j,')=NTERMS_DENS (=',ncutoffDens_overall,')'
                          ncutoffDens(i,j)=ncutoffDens_overall
                       elseif (ncutoffs.gt.0) then
                          !...or if NTERMS_DENS is not specified use the
                          !'master' NTERMS value.
                          write(*,'(A13,I1,A1,I1,A11,I2,A1)') &
                ' NTERMS_DENS(',i,';',j,')=NTERMS (=',ncutoffs,')'
                          ncutoffDens(i,j)=ncutoffs
      !                 else
      !                    write(*,'(A16,I1,A1,I1,A48)') &
      !       ' No NTERMS_DENS(',i,';',j,') in settings file, MUST be specified, STOPPING.'
      !                    stop
                       endif
                    endif
                 endif

                 if ((readParasOnebyone.eqv..false.).and.(startparas.eq.2)) then
                    !Update the freep array so code optimizes the relevant
                    !electron density parameters.
                    do l=0,lmax
                       do ipara=paraCnt,paraCnt+2*ncutoffDens(i,j)-1
                          freep(ipara)=2
                       enddo
                       paraCnt=paraCnt+12
                    enddo
                 endif
                        
              else
                 paraCnt=paraCnt+12*lm1
              endif
           enddo
        enddo 
        
        !Read-in the number of terms, per species, to be used in the embedding
        !function
        searchfor="NTERMS_EMB"
        call getsettingspara(searchfor,found,foundvariable)
        if (foundvariable) then
           read(found,*) nembterms
           if (readParasOnebyone) then
             print *,'ERROR: NTERMS_EMB tag conflicts with PARASTOOPT. Please'
             print *,'remove one and restart, stopping.'
             stop
           endif
           write(*,'(A12,I1,A21)') ' NTERMS_EMB=',nembterms,' (from settings file)'
           if (nembterms.eq.0) then
              print *,' ERROR: NTERMS_EMB=0, stopping.'
              stop
           endif
           if (nembterms.gt.4) then
              print *,' ERROR: NTERMS_EMB>4, stopping.'
              stop
           endif
           startparas=2
       else
          nembterms=0
!          if (readParasOnebyone.eqv..false.) then
!             write(*,'(A60)') &
!            ' No NTERMS_EMB in settings file, MUST be specified, STOPPING.'
!             stop
!          endif
        endif

        if (readParasOnebyone.eqv..false.) then
           !Based on the NTERMS_EMB setting, fill in relevant parts of freep array 
           !to tell code to optimize the relevant embedding function parameters.
           if (nembterms.ge.1) then
              freep(2*m3+lm1*m1+12*lm1*m2+1:2*m3 &
         +lm1*m1+12*lm1*m2+m1)=2
           endif
           if (nembterms.ge.2) then
              freep(2*m3+lm1*m1+12*lm1*m2+1+m1:2*m3 &
         +lm1*m1+12*lm1*m2+2*m1)=2
           endif
           if (nembterms.ge.3) then
              freep(2*m3+lm1*m1+12*lm1*m2+1+2*m1:2*m3 &
         +lm1*m1+12*lm1*m2+3*m1)=2
           endif
           if (nembterms.eq.4) then
              freep(2*m3+lm1*m1+12*lm1*m2+1+3*m1:2*m3 &
         +lm1*m1+12*lm1*m2+4*m1)=2
           endif
        endif

        !As described above, allows number of terms in pair-potentials to be
        !specified independent of number of terms in electron-densities
        nterms_pp_overall=0
        searchfor="NTERMS_PP"
        call getsettingspara(searchfor,found,foundvariable)
        if (foundvariable) then
           if (readParasOnebyone) then
             print *,'ERROR: NTERMS_PP tag conflicts with PARASTOOPT. Please'
             print *,'remove one and restart, stopping.'
             stop
           endif
           read(found,*) nterms_pp_overall
           if (verbose) then
              write(*,'(A11,I1,A21)') ' NTERMS_PP=',nterms_pp_overall,' (from settings file)'
           endif
           startparas=2
        endif

        paraCnt=2*m3+(4+lm1)*m1+12*lm1*m2+1
        do i=1,maxspecies
           do j=1,maxspecies
              if (j.ge.i) then
                 !Check if NTERMS_PP(1,2) etc are supplied, which as for
                 !electron densities allows for number of terms to be different
                 !for different pairs of species.
                 write(string2,'(I1)') i
                 write(string3,'(I1)') j
                 searchfor="NTERMS_PP("//string2//";"//string3//")"
                 call getsettingspara(searchfor,found,foundvariable)
                 if (foundvariable) then
                    if (readParasOnebyone) then
                       write(*,'(A8,A14,A31)') ' ERROR: ',trim(searchfor),' tag conflicts with PARASTOOPT.'
                       print *,'Please remove one and restart, stopping.'
                       stop
                    endif
                    read(found,*) nterms_pp(i,j)
                    write(*,'(A11,I1,A1,I1,A2,I2,A21)') &
          ' NTERMS_PP(',i,';',j,')=',nterms_pp(i,j),' (from settings file)'
                    startparas=2
                 else
                    if ((readParasOnebyone.eqv..false.).and.(startparas.eq.2)) then
                       if (nterms_pp_overall.gt.0) then
                          write(*,'(A19,I1,A1,I1,A16,I2,A1)') &
                ' setting NTERMS_PP(',i,';',j,')=NTERMS_PP (=',nterms_pp_overall,')'
                          nterms_pp(i,j)=nterms_pp_overall
                       elseif (ncutoffs.gt.0) then
                          write(*,'(A19,I1,A1,I1,A11,I2,A1)') &
                ' setting NTERMS_PP(',i,';',j,')=NTERMS (=',ncutoffs,')'
                          nterms_pp(i,j)=ncutoffs
                       else
                          nterms_pp(i,j)=0
           !              write(*,'(A14,I1,A1,I1,A48)') &
           ! ' No NTERMS_PP(',i,';',j,') in settings file, MUST be specified, STOPPING.'
           !              stop
                       endif
                    endif
                 endif

                 !Fill freep array to tell MEAMfit which pair-potential parameters to
                 !optimize
                 if ((readParasOnebyone.eqv..false.).and.(startparas.eq.2)) then
                    if (nterms_pp(i,j).gt.0) then
                       do ipara=paraCnt,paraCnt+2*nterms_pp(i,j)-1
                          freep(ipara)=2
                       enddo
                    endif
                 endif

              endif
              paraCnt=paraCnt+32
           enddo
        enddo

        !Whether to optimize the energy constant (per species) - see paper. This
        !should normally be true, except for example where another potential is
        !used as a starting point for optimization of, eg, a cross-potential.
        searchfor="OPT_ENCONST"
        call getsettingspara(searchfor,found,foundvariable)
        if (foundvariable) then
           if (readParasOnebyone) then
             print *,'ERROR: OPT_ENCONST tag conflicts with PARASTOOPT. Please'
             print *,'remove one and restart, stopping.'
             stop
           endif
           read(found,*) opt_enconst
           write(*,'(A13,L1,A21)') ' OPT_ENCONST=',opt_enconst,' (from settings file)'
           startparas=2
        else
           if (readParasOnebyone.eqv..false.) then
              if (verbose) then
                 print *,'No OPT_ENCONST in settings file, using default (=.false.)'
              endif
              opt_enconst=.false.
           endif
        endif

        !Whether to optimize the 'meam t' coefficients appearing in the
        !background density (total electron density)
        searchfor="OPT_MEAMTAU"
        call getsettingspara(searchfor,found,foundvariable)
        if (foundvariable) then
           if (readParasOnebyone) then
             print *,'ERROR: OPT_ENCONST tag conflicts with PARASTOOPT. Please'
             print *,'remove one and restart, stopping.'
             stop
           endif
           read(found,*) opt_meamtau
           if (verbose) then
              write(*,'(A12,L1,A21)') ' OPT_MEAMTAU',opt_meamtau,' (from settings file)'
           endif
        else
           if (readParasOnebyone.eqv..false.) then
              if (startparas.eq.2) then
                 if (lmax.eq.0) then
                    if (verbose) then
                       print *,'No OPT_MEAMTAU in settings file, using default for TYPE=EAM (=.false.)'
                    endif
                    opt_meamtau=.false.
                 elseif (lmax.eq.3) then
                    print *,'No OPT_MEAMTAU in settings file, using default for TYPE=MEAM (=.true.)'
                    opt_meamtau=.true.
                 endif
              else
                 if (verbose) then
                    print *,'No OPT_MEAMTAU in settings file, using default (=.false., since no other optimizable'
                    print *,'parameters specified)'
                 endif
                 opt_meamtau=.false.
              endif
           endif
        endif

        if (readParasOnebyone.eqv..false.) then

           if (opt_meamtau) then
              do ipara=2*m3+m1+1,2*m3+lm1*m1
                 freep(ipara)=2
              enddo
           endif

           if (opt_enconst) then
              !Only one enconst varied by default (assumes same number of
              !species per configuration) Andy: this will need to be changed for
              !the case where there are more than two species.
              freep(m4+2*m2+1)=2
            ! do ipara=m4+2*m2+1,m4+2*m2+m1
            !     freep(ipara)=2
            ! enddo
           endif

        endif

        deallocate(ncutoffDens,ncutoffPairpot,nterms_pp)

       endif

      !if ((startparas.eq.1).and.(readpotfile.eqv..false.)) then
      !   print *,'ERROR: Potential parameter read-in specified (NOOPT=true or'
      !   print *,'-noopt/-no flag at command-line), yet no potential filename'
      !   print *,'supplied, STOPPING.'
      !   stop
      !endif

    endif

    !Tally the total number of parameters to be randomly generated. Andy: this
    !does not appear to be used at the moment, consider excising.
    if ((nsteps.gt.1).or.(startparas.eq.2)) then
       nfreep=0
       nRandGen=0
       do i=1,np
           if ((freep(i).eq.1).or.(freep(i).eq.2)) then
              nfreep=nfreep+1
           endif
           if (freep(i).eq.2) then
               if (startparas.eq.2) then
                   nRandGen=nRandGen+1
               else
                   !setting startparas=1 overrides the other optimization settings
                   freep(i)=1
               endif
           endif
           if (freep(i).eq.3) then
               if (startparas.eq.2) then
                   nRandGen=nRandGen+1
               else
                   !setting startparas=1 overrides the other optimization settings
                   freep(i)=0
               endif
           endif

        !  if (freep(i).ge.1) then
        !      nfreep=nfreep+1
        !      if (freep(i).eq.2) then
        !          if (startparas.eq.2) then
        !              nRandGen=nRandGen+1
        !          else
        !              !setting startparas=1 overrides the other optimization settings
        !              freep(i)=1
        !          endif
        !      endif
        !  endif

       enddo
    endif

! code fragment (from earlier); didn't work because I couldn't find a way
! to put the value into the varable with the name 'NamePara(i)'.

     ! character*10, NamePara(2)
     ! character*1, TypePara(2)
     ! character*10, DfltPara(2)
     ! logical ReqPara(2)

     ! data NamePara/'VERBOSE','SEED'/
     ! data TypePara/'L','I'/
     ! data DfltPara/'.TRUE.','0'/
     ! data ReqPara/'F','F'/

     ! nPara=2

     ! verbose=.false.
     ! do i=1,nPara
     !    searchfor=trim(NamePara(i))
     !    call getsettingspara(searchfor,found,typePara(i),foundvariable)
     !    !print *,'looked for:',searchfor,' of type ',typePara(i)
     !    !print *,'found:',foundvariable
     !    if (foundvariable.eqv..true.) then
     !       read(found,*) NamePara(i)
     !       print *,'found ',searchfor,'=',NamePara(i),' in settings file.'
     !       print *,'should have been assigned to verbose:',verbose
     !       stop
     !    else
     !       if (ReqPara(i).eqv..true.) then
     !          print *,'ERROR: '
     !       else
     !       endif
     !    endif
     ! enddo
     ! stop

end subroutine readsettings
subroutine resetxyz

    !Reset xyz files

    use m_filenames
    use m_geometry

    character*80 poscarheader,poscarspecies
    integer i,j
    real(8) latticeconst,poscarcell(9)

    print *,'check that the following are the same as originals'
    do i=1,nstruct
        open(unit=1,file=trim(poscarfiles(i)))
        read(1,'(a)') poscarheader
        read(1,*) latticeconst
        read(1,*) poscarcell
        read(1,'(a)') poscarspecies
        close(1)
        rewind(unit=1)
        write(1,'(a)') trim(poscarheader)
        write(1,*) latticeconst
        write(1,'(3f12.6)') poscarcell
        write(1,'(a)') trim(poscarspecies)
        write(1,'(a)') 'Cartesian'
        do j=1,gn_forces(i) !or gn_inequivalentsites(i) ?
            write(1,*) gxyz_backup(1:3,j,i)/latticeconst
        enddo
        close(1)
    enddo
    print *,'stopping after resetxyzx...'
    stop
end subroutine resetxyz
subroutine rmforcenoise(en1,en2,dist,force,forcenew)

    !----------------------------------------------------------------------
    !
    !     Determine to how many s.f. the force is accurate to (assuming
    !     force calculation using finite difference method: force =
    !     - ( en2 - en1 ) / dist), and then truncate the remaining s.f.'s
    !
    !     Andrew Duff, 2014
    !
    !----------------------------------------------------------------------

    use m_generalinfo

    implicit none

    integer orderendiff,orderen,orderfrc,frcprc,numdptomv,nint,tmp
    real(8) en1,en2,dist,force,forcenew

        ! write(*,'(A18,F30.20)') 'original energy: ',en1
        ! write(*,'(A18,F30.20)') 'displaced energy:',en2
        ! write(*,'(A6,F30.20,A8,F30.20)') 'diff=',en2-en1,', dist=',dist
    orderendiff=floor(log10(abs(en2-en1))) !log10 gets the order,
    !floor rounds down
        !  write(*,*) '...of order: O(',orderendiff,')'
    orderen=floor(log10(abs(en1)))
        !  write(*,*) 'energy of order: O(',orderen,')'
        !  write(*,*) 'therefore difference in energy occurs at ', &
        !         -orderendiff+orderen+1,'th s.f. of energy'
        !  write(*,*) 'nsigdig for all calcs=',nsigdig
    if (-orderendiff+orderen+1.gt.nsigdig) then
          !       write(*,*) 'diff in en thus occurs greater than available', &
          !              ' precision, thus settings force=0'
        forcenew=0d0
    else
        frcprc=nsigdig+orderendiff-orderen
         !        write(*,*) 'diff, and thus force, therefore only precise to ', &
         !               frcprc
        !the following moves all s.f. to left of d.p., converts to
        !integer, then moves the digits back to their original
        !positions
        !(thus settings all those digits > no. s.f. to zero)
        orderfrc=floor(log10(abs(force)))
        numdptomv=frcprc-orderfrc-1
        !         write(61,*) 'moving digits in force left ',numdptomv,' digits'
        !         forcenew=floor(force*dble(10.d0**numdptomv))/
        !     +            dble(10.d0**numdptomv)
        tmp=nint(force*dble(10.d0**numdptomv))
        forcenew=dble(tmp)/dble(10.d0**numdptomv)
    endif
    ! write(*,'(A18,F30.20)') 'original force:  ',force
    ! write(*,'(A18,F30.20)') 'original force:  ',forcenew
    ! stop
end subroutine rmforcenoise
subroutine screening_ij

    !---------------------------------------------------------------c
    !
    !     Fills the module m_screening with the S_ij parameters
    !     (NOTE: the screening parameters Cmin and Cmax are species
    !     dependent)
    !
    !     Called by:     meamenergy
    !     Calls:         distance2
    !     Arguments:     istr,gn_inequivalentsites,gn_neighbors,
    !                    cmin_cmax_zero,gn_neighbors,gspecies,gneighborlist,
    !                    cmax,cmin
    !     Returns:       screening
    !     Files read:    -
    !     Files written: -
    !
    !     Marcel Sluiter, Feb 8 2006
    !
    !---------------------------------------------------------------c

    use m_atomproperties
    use m_geometry
    use m_meamparameters
    use m_screening   !S_ij

    implicit none

    integer i,k,j,jj,kk,nni,iaux(1),isp,jsp,ksp
    real(8) s_ij,alpha_ik,alpha_jk,aux,c_ikj,y_ikj,y,s_ikj

    if(allocated(screening)) deallocate(screening)
    iaux=maxval(gn_neighbors(1:gn_inequivalentsites(istr),istr))
    allocate(screening(iaux(1),gn_inequivalentsites(istr)))
    if (cmin_cmax_zero) then

        do i=1,gn_inequivalentsites(istr)

            nni=gn_neighbors(i,istr)
            do jj=1,nni
                screening(jj,i)=1.d0
            enddo
        enddo

    else

        do i=1,gn_inequivalentsites(istr)
            isp=gspecies(i,istr)
            nni=gn_neighbors(i,istr)
            do jj=1,nni
                j=gneighborlist(jj,i,istr)
                jsp=gspecies(j,istr)
                s_ij=1d0
                kk=1
                do while(s_ij.ne.0d0.and.kk.le.nni)
                    if(kk.ne.jj) then !avoid k=j
                        k=gneighborlist(kk,i,istr)
                        ksp=gspecies(k,istr)
                        if ((cmin(isp,ksp,jsp).ne.0d0).and. &
                            (cmax(isp,ksp,jsp).ne.0d0)) then
                            alpha_ik=distance2(i,k)/distance2(i,j)
                            alpha_jk=distance2(j,k)/distance2(i,j)
                            print *,distance2(i,k),dist2str(kk,i,0,0,istr)
                            aux=(alpha_ik-alpha_jk)**2
                            c_ikj=(2d0*(alpha_ik+alpha_jk)-aux-1d0)/ &
                                (1d0-aux)
                            if((1.d0-aux).lt.0.d0) then
                                s_ikj=1d0
                            else
                                y_ikj=(c_ikj-cmin(isp,ksp,jsp))/ &
                                    (cmax(isp,ksp,jsp)-cmin(isp,ksp,jsp))
                                y=min(1d0,y_ikj) !for y>1 result same as for
                                !y=1
                                s_ikj=0d0
                                if(y.gt.0d0) s_ikj=(1d0-(1d0-y)**4)**2
                            endif
                            s_ij=s_ij*s_ikj
                            !                        if (s_ikj.ne.1d0) then
                            !                           print *,
                            !     +     'istr,i,distance2(i,k),distance2(j,k),',
                            !     +     'distance2(i,j),alpha_ik,alpha_jk,c_ikj,cmin(isp,ksp,jsp),',
                            !     +        'cmax(isp,ksp,jsp),y_ikj,s_ikj,isp,jsp,ksp'
                            !                          print *,istr,i,distance2(i,k),distance2(j,k),
                            !     +        distance2(i,j),alpha_ik,alpha_jk,c_ikj,cmin(isp,ksp,jsp),
                            !     +        cmax(isp,ksp,jsp),y_ikj,s_ikj,isp,jsp,ksp
                            !c                           stop
                            !                        endif
                        endif
                    endif
                    kk=kk+1
                enddo
                screening(jj,i)=s_ij
            enddo

        enddo

    endif

end subroutine screening_ij
subroutine setupdefaultlmax

    !--------------------------------------------------------------c
    !
    !     Use values of 'maxspecies' and 'lmax' to set up the    
    !     variables to be used to reference the potential parameter
    !     array, p().
    !
    !     Called by:     program MEAMfit
    !     Calls:         
    !     Arguments:     
    !     Returns:       
    !     Files read:    
    !     Files written: -
    !
    !     Andrew Duff 2014
    !
    !--------------------------------------------------------------c

    use m_meamparameters

    implicit none

    integer, parameter :: lmaxDefault = 0 ! Number of angular momentum terms to set
                                          ! up for in the electron density terms

    lmax=lmaxDefault

end subroutine setupdefaultlmax
subroutine setupfilenames

    !-------------------------------------------------------------c
    !
    !     Assigns filenames to the filename variables. Also sets
    !     up the array 'strucnames', which contains abbreviated
    !     forms of the vasprun.xml/POSCAR filenames, specifically, 
    !     only the parts of the filenames pertinant to defining the
    !     structure. This is explained by way of example. Say we
    !     have a filename vasprun_Coct_pot_relaxed.xml, then the 
    !     part of this filename which goes into the relevant
    !     'strucnames' array element is 'Coct'. The vasprun_ is
    !     always removed, and if either 'pot', 'EAM' or 'MEAM'
    !     are identified, then these (and everything to the right
    !     of them) are also removed.
    !
    !     Called by:     program MEAMfit
    !     Calls:         createfitdbse, readFitdbseLine,
    !                 checkifnumber, getnumconfigs 
    !     Returns:       nposcarfiles,nstruct,startparameterfile,
    !                 crystalstrucfile,settingsfile,
    !                 fitdbse, poscarfiles,outcarfiles,strucnames
    !     Files read:    fitdbse
    !
    !     Andy Duff, Feb 2008
    !
    !-------------------------------------------------------------c

    use m_filenames
    use m_geometry
    use m_datapoints
    use m_optimization

    implicit none

    logical num,ensureminimum_tocopy,freeenergy_tocopy,vasprun_tocopy,exist
    integer i,j,ifile,istruc,nlimits,optimizeforce_tocopy
    integer limits(10),step(10)
    character*80 line
    character*80 tmp,tmp2,string,string2,string3,string4
    character*80 poscarfiles_tocopy,outcarfiles_tocopy
    real(8) weights_tocopy,double

    crystalstrucfile='crystal_struc.xyz'
    settingsfile='settings'
    fitdbse='fitdbse'
    lookuptablefile='lookuptablefile'
 
    !Check if there is a fitdbse file
    inquire(file=trim(fitdbse),exist=exist)
    if (exist.eqv..false.) then
       call createfitdbse
    endif
 
    !Read in list of poscar files
    open(unit=1,file=trim(fitdbse),status='old')
    read(1,*) nposcarfiles
    if (nposcarfiles.lt.1) then
        print *,'Error: no. of poscarfiles must be greater than'
        print *,'or equal to 1'
        stop
    endif
    !Scan through fitdbse once to see how many ionic structures we have to fit
    !to (this will not necessarily equal the number of VASP files)
    nstruct=0
    do ifile=1,nposcarfiles
 
        read(1,'(A80)') line
        call readFitdbseLine(line,tmp,tmp2,string,string4,double)
        call checkifnumber(string,num)
        if (num.eqv..true.) then
           call getnumconfigs(string,limits,nlimits,step)
        !  print *,'string=',string
        !  print *,'limits=',limits
        !  print *,'nlimits=',nlimits
        !  print *,'step=',step
        !  stop
           do i=1,nlimits,2
              do j=limits(i),limits(i+1),step(i)
                 nstruct=nstruct+1
              enddo
           enddo
        else
           nstruct=nstruct+1
        endif
    enddo

    rewind(1)
    read(1,*)
    allocate(poscarfiles(nstruct),outcarfiles(nstruct), &
        strucnames(nstruct),optimizeforce(nstruct), &
        ensureminimum(nstruct),weights(nstruct), &
        rlxstruc(nstruct),optimizeforce_backup(nstruct), &
        freeenergy(nstruct),nconfig(nstruct), &
        vasprun(nstruct))
    istruc=1
  
    print *
    print *,'Initializing Fitting database'
    print *,'-----------------------------'
    print *,'File | Configs to fit | Quantity to fit | Weights'

    forcefit=.false.  
    do ifile=1,nposcarfiles
  
       !New parsing code to improve transferability:
       read(1,'(A80)') line
       call readFitdbseLine(line,poscarfiles(istruc),outcarfiles(istruc),string,string4,weights(istruc))

       string2=poscarfiles(istruc)
       if (string2(1:7).eq."vasprun") then
          vasprun(istruc)=.true.
          outcarfiles(istruc)=poscarfiles(istruc)
          write(string3,'(i10)') nconfig(istruc)
          strucnames(istruc)=trim(poscarfiles(istruc))//trim(string3)
       else
          print *,'ERROR: vasprun files in fitdbse file must start with the word vasprun.'
          print *,'Stopping.'
          stop
          !The following is a relic from when POSCAR/OUTCAR combinations were allowed
          vasprun(istruc)=.false.
       endif

       !second columns: n=no ionic relaxation and read in one ionic
       !config from file; r=relax structure using potential; numerical
       !sequence (e.g, 1 or 1,2 or 1-5 or 1-5,6-10): read in those
       !ionic configurations from the file.
       call checkifnumber(string(1:1),num)
       if (num.eqv..false.) then
          nconfig(istruc)=1
          !VASP file has only one ionic config to be read in. Should the ionic
          !positions be relaxed (according to the potential) during optimization?
          if (string(1:1).eq.'N'.or.string(1:1).eq.'n') then
              rlxstruc(istruc)=0 !Do not relax structure using potential
          elseif (string(1:1).eq.'R'.or.string(1:1).eq.'r') then
              rlxstruc(istruc)=1 !Relax structure using potential
          else
              write(*,'(A26,A1,A37)') 'ERROR: Unexpected letter: ',string(1:1), &
                              ', in 3rd column of fitdbse. STOPPING.'
              stop
          endif
       endif
  
       !If the third column contains an 'F' of 'f' then set
       !optimizeforce=.true. for this file, which signals that the code
       !should optimize w.r.t the forces from this outcar file;
       !otherwise set optimizeforce=.false. and therefore optimize
       !w.r.t the energy from this file.
       if (string4(1:2).eq.'fo'.or.string4(1:2).eq.'Fo'.or. &
           string4(1:2).eq.'fO'.or.string4(1:2).eq.'FO') then
           optimizeforce(istruc)=1
           forcefit=.true.
      !elseif(snglechar2(1:1).eq.'G'.or.snglechar2(1:1).eq.'g') then
      !    optimizeforce(istruc)=2
      !    forcefit=.true.
       else
           optimizeforce(istruc)=0
          !if(snglechar2(1:1).eq.'M'.or.snglechar2(1:1).eq.'m') then
          !    ensureminimum(istruc)=.true.
          !else
           ensureminimum(istruc)=.false.
          !endif
           if (string4(1:2).eq.'fr'.or.string4(1:2).eq.'Fr'.or. &
               string4(1:2).eq.'fR'.or.string4(1:2).eq.'FR') then
               freeenergy(istruc)=.true.
           elseif (string4(1:2).eq.'e0'.or.string4(1:2).eq.'E0') then
               freeenergy(istruc)=.false.
           else
               print *,'ERROR: Incorrect entry in 3rd column of fitdbse file.'
               print *,'Expecting E0, FR or FO. Stopping'
               stop
           endif
       endif

       if (optimizeforce(istruc).gt.0) then
          string2='Force'
       else
          if (freeenergy(istruc)) then
             string2='Free energy'
          else
             string2='Energy (E0)'
          endif
       endif
       if (vasprun(istruc).eqv..false.) then
          print *,trim(poscarfiles(istruc)),' ',trim(outcarfiles(istruc)), &
                  ' ',trim(string),' ',trim(string2),' ',weights(istruc)
       else
          print *,trim(poscarfiles(istruc)),' ',trim(string),' ', &
                  trim(string2),' ',weights(istruc)
       endif


       if (num.eqv..true.) then
          rlxstruc(istruc)=0
          vasprun_tocopy=vasprun(istruc)
          poscarfiles_tocopy=poscarfiles(istruc)
          outcarfiles_tocopy=outcarfiles(istruc)
          optimizeforce_tocopy=optimizeforce(istruc)
          ensureminimum_tocopy=ensureminimum(istruc)
          freeenergy_tocopy=freeenergy(istruc)
          weights_tocopy=weights(istruc)
          call getnumconfigs(string,limits,nlimits,step)
          do i=1,nlimits,2
             !print *,'i=',i,' going from: ',limits(i),'to ',limits(i+1)
             do j=limits(i),limits(i+1),step(i)
                vasprun(istruc)=vasprun_tocopy
                poscarfiles(istruc)=poscarfiles_tocopy
                outcarfiles(istruc)=outcarfiles_tocopy
                rlxstruc(istruc)=0
                optimizeforce(istruc)=optimizeforce_tocopy
                ensureminimum(istruc)=ensureminimum_tocopy
                freeenergy(istruc)=freeenergy_tocopy
                weights(istruc)=weights_tocopy
                nconfig(istruc)=j
                if (nconfig(istruc).lt.10) then
                   write(string3,'(I1)') nconfig(istruc)
                elseif (nconfig(istruc).lt.100) then
                   write(string3,'(I2)') nconfig(istruc)
                elseif (nconfig(istruc).lt.1000) then
                   write(string3,'(I3)') nconfig(istruc)
                elseif (nconfig(istruc).lt.10000) then
                   write(string3,'(I4)') nconfig(istruc)
                endif
                strucnames(istruc)=trim(poscarfiles_tocopy)//"_"//trim(string3)
                istruc=istruc+1
             enddo
          enddo
       else
          istruc=istruc+1
       endif
  
     enddo
     write(*,'(A13,I5,A40)') ' (fitting to ',nstruct,' atomic configurations across all files)'
     close(unit=1)

end subroutine setupfilenames
subroutine setupnntables

    !----------------------------------------------------------------------c
    !
    !     Set-up interatomic seperation arrays, to avoid having to
    !     recalculate the interatomic seperations each time they are used.
    !     Note: array are also initialized for the force calculations. In
    !     these cases, seperations are also recorded for the configurations
    !     where each atom, taken in turn, is displaced in the 3 cartesian
    !     directions.
    !
    !     Called by:     program MEAMfit
    !
    !     Andy Duff, Oct 2013
    !
    !----------------------------------------------------------------------c

    use m_geometry
    use m_datapoints

    implicit none

    integer i,j,jj,nni,idisp,cart
    real(8), parameter:: dist=1d-8 !The distance moved in the x,y and
    !z directions for the force-calculation

    !Store all interatomic distances in arrays to speed up code (also squares
    !and cubes of distances, and cartesian components of seperations)
    do istr=1,nstruct
        do i=1,gn_inequivalentsites(istr)
            nni=gn_neighbors(i,istr)
            do jj=1,nni
                j=gneighborlist(jj,i,istr)
                diststr(jj,i,0,0,istr)=distance(i,j)
                dist2str(jj,i,0,0,istr)=distance2(i,j)
                dist3str(jj,i,0,0,istr)= &
                    diststr(jj,i,0,0,istr)*dist2str(jj,i,0,0,istr)
                dxstr(jj,i,0,0,istr)=gxyz(1,j,istr)-gxyz(1,i,istr)
                dystr(jj,i,0,0,istr)=gxyz(2,j,istr)-gxyz(2,i,istr)
                dzstr(jj,i,0,0,istr)=gxyz(3,j,istr)-gxyz(3,i,istr)
            enddo
        enddo

        !Also store seperations for force calculations, I.e., the seperations
        !between atoms in the displaced structures (there will be 3*natoms of
        !these, where natoms is the number of atoms in the structure)
        if ((optimizeforce(istr).gt.0).and.(fastForce)) then
            !            print *,'setupnntables, istr=',istr
            do idisp=1,gn_forces(istr) !Atom to displace
                if (optforce(idisp,istr).eqv..true.) then
                    !                 print *,'idisp=',idisp,', istr=',istr
                    do cart=0,3 !loop over cartesian coords
                        if (cart.gt.0) then
                            gxyz(cart,idisp,istr)=gxyz(cart,idisp,istr)+dist
                        endif
                        do i=1,gn_inequivalentsites(istr)
                            nni=gn_neighbors(i,istr)
                            do jj=1,nni
                                j=gneighborlist(jj,i,istr)
                                diststr(jj,i,idisp,cart,istr)=distance(i,j)
                                dist2str(jj,i,idisp,cart,istr)= &
                                    distance2(i,j)
                                dist3str(jj,i,idisp,cart,istr)= &
                                    diststr(jj,i,idisp,cart,istr)*dist2str(jj,i,idisp,cart,istr)
                                dxstr(jj,i,idisp,cart,istr)= &
                                    gxyz(1,j,istr)-gxyz(1,i,istr)
                                dystr(jj,i,idisp,cart,istr)= &
                                    gxyz(2,j,istr)-gxyz(2,i,istr)
                                dzstr(jj,i,idisp,cart,istr)= &
                                    gxyz(3,j,istr)-gxyz(3,i,istr)
                            enddo
                        enddo
                        if (cart.gt.0) then
                            gxyz(cart,idisp,istr)=gxyz(cart,idisp,istr)-dist
                        endif
                    enddo
                endif
            enddo
        endif

    enddo

end subroutine setupnntables
subroutine setuppositivep

    !----------------------------------------------------------------------c
    !
    !     Set up positivep array, which determines which parameters in the
    !     p() array are to remain positive.
    !
    !     Called by: optimizeparameters
    !     Calls: -
    !     Returns: positivep
    !     Files read: -
    !     Files written: -
    !
    !     Andrew Duff, 2014
    !
    !----------------------------------------------------------------------c


    use m_generalinfo
    use m_meamparameters
    use m_optimization
    use m_atomproperties

    implicit none

    integer i,j,k,paracounter,isp,jsp,spcCnt

    positivep=.false.

    !Atomic densities: cut-off radii should be positive
    paracounter=2*m3+lm1*m1
    do i=1,maxspecies
        do j=1,maxspecies
            if ((i.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                if (typethi(i,j).eq.1) then
                   do k=paracounter+2,paracounter+12*lm1,2
                      positivep(k)=.true.
                   enddo
                else
                   print *,'typethi(',i,',',j,') (=',typethi(i,j),') not supported, STOPPING.'
                   stop
                endif
            endif
            paracounter=paracounter+12*lm1
        enddo
    enddo

    !Embedding function
    !------------------
    if (embfunctype.eq.1) then
       !The rho^1/2 coefficient should always be positive (and finally
       !multiplied by a minus sign)
       positivep(2*m3+lm1*m1+12*lm1*m2+1: &
                 2*m3+lm1*m1+12*lm1*m2+m1)=.true.
    else
       print *,'embfunctype=2 not supported, stopping.'
       stop
    endif

    !Pairpotential
    !-------------
    !Cut-offs should be positive
    spcCnt=0
    do isp=1,m1
        do jsp=1,m1
            if (jsp.ge.isp) then
               if (typepairpot(isp,jsp).eq.2) then
                  do i=2*m3+(4+lm1)*m1+12*lm1*m2+2,2*m3+(4+lm1)*m1+12*lm1*m2+32,2
                      positivep(i+32*spcCnt)=.true.
                  enddo
               endif
            endif
            spcCnt=spcCnt+1
        enddo
    enddo

end subroutine setuppositivep
subroutine setuprandomseed(iseed)

    !-------------------------------------------------------------------c
    !
    !     Initializes the random seed, either based on the current
    !     time or a user supplied seed from the settings file.
    !
    !     Andrew Duff 2015
    !
    !-------------------------------------------------------------------c
  

    implicit none

    integer seed(12) !Changed from 8
    character*8 charac1
    character*10 charac2
    character*5 charac3
    integer values(8),iseed,i
    real(8) tmp

    if (iseed.eq.0) then
       !Initialize random seed using data and time
       call date_and_time(charac1, charac2, charac3, values)
       do i=1,8
          seed(i)=values(8)*values(i)
       enddo
    else
       !Initialize random seed using integer provided in settings file
       do i=1,8
          seed(i)=iseed*dble(i)
       enddo
    endif
   
   call random_seed(put=seed)

   !test
   !call RANDOM_NUMBER(tmp)
   !print *,tmp
   !call RANDOM_NUMBER(tmp)
   !print *,tmp
   !call RANDOM_NUMBER(tmp)
   !print *,tmp
   !call RANDOM_NUMBER(tmp)
   !print *,tmp
   !call RANDOM_NUMBER(tmp)
   !print *,tmp
   !stop
   !!:::::::

end subroutine setuprandomseed
subroutine setupsplines

    use m_meamparameters
    use m_optimization
    use m_generalinfo

    implicit none

    integer i
    real(8) pairpotAcklandFeFe

    !Setup spline for V_XX
    do i=1,narrpairpotXX
        pairpotXXarr(i)=pairpotAcklandFeFe(r_pairpotXXarr(i))
    enddo
    call spline(r_pairpotXXarr,pairpotXXarr,narrpairpotXX, &
        0d0,0d0,secderpairpotXX)

    !Setup spline for thi_XX
    do i=1,narrthiXX
        call raddens(r_thiXXarr(i),1,1,thiXXarr(i))
    enddo
    call spline(r_thiXXarr,thiXXarr,narrthiXX, &
        0d0,0d0,secderthiXX)

end subroutine setupsplines
subroutine spline(X,Y,n,YP1,YPN,y2)


    !-----------------------------------------------------------------
    !     From Numerical Recipes
    !     DESCRIPTION: fits a cubic spline to a supplied function
    !                  y(X).
    !     INPUT: A function y(X) stored in the array
    !            y(1:N_ARR) defined at points X(1:n).
    !            Points from 1 to
    !            N are included. dY/dX must also be specified at
    !            X(1) and X(N) and are passed to the subroutine as
    !            YP1 and YPN.
    !     OUTPUT: Subroutine calculates d^2Y/dx^2 at all x points
    !             and returns it in the array Y2(1:N)
    !     INSTRUCTIONS: run this subroutine once to calculate Y2,
    !                   then 'splint' can be run to calculate the
    !                   value of the interpolated Y(X) at any X. Or
    !                   'spline_integrate' can be run to integrate
    !                   the spline.
    !                   to fit a natural spline (d^2Y/dx^2=0 at
    !                   X(1) ), set YP1 or YPN to > 10^30
    !-----------------------------------------------------------------

    implicit none
    integer i,n,nmax
    real *8 YP1,YPN,X(n),Y(n),y2(n)
    parameter (nmax=120000)!3000
    INTEGER k
    REAL *8 p,qn,sig,un,u(nmax)
    !      print *,'n=',n,', nmax=',nmax
    if (n.gt.nmax) then
        print *,'nmax in spline needs increasing. Stopping.'
        stop
    endif
    !      if (YP1.gt..99e30) then    !The lower boundary condition is set
    !                                 !either to be
    y2(1)=0.d0
    u(1)=0.d0
    !      else                      !or else to have a specified first
    !         y2(1)=-0.5d0           !derivative.
    !         u(1)=(3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    !      endif
    do i=2,n-1             !This is the decomposition loop of the
        sig=(X(i)-X(i-1))/(X(i+1)-X(i-1)) !tridiagonal algorithm. y2
        p=sig*y2(i-1)+2.d0  !and u are used fo temporary storage of
        y2(i)=(sig-1.d0)/p  !the decomposed factors.
        u(i)=(6.d0*((Y(i+1)-Y(i))/(X(i+1)-X(i))-(Y(i)-Y(i-1)) &
            /(X(i)-X(i-1)))/(X(i+1)-X(i-1))-sig*u(i-1))/p
    enddo
    !      if (YPN.gt..99e30) then    !The upper boundary condition is set
    qn=0.d0                  !either to be
    un=0.d0
    !      else                      !or else to have a specified first
    !         qn=0.5d0               !derivative.
    !         un=(3.d0/(X(n)-X(n-1)))*(YPN-(Y(n)-Y(n-1))/
    !     +        (X(n)-X(n-1)))
    !      endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
    do k=n-1,1,-1          !This is the backsubstitution loop of the
        y2(k)=y2(k)*y2(k+1)+u(k) !tridiagonal algorithm.
    enddo

end subroutine spline
subroutine splint(XA,YA,Y2A,N,N_ARR,index,X,Y)


    !----------------------------------------------------------------c
    !     From Numerical Recipes                                     c
    !     DESCRIPTION: given arrays YA(1:N) and XA(1:N) which        c
    !                  tabulate a function and Y2A(1:N) as           c
    !                  calculated by 'spline', the subroutine        c
    !                  calculates the value of the spline at X and   c
    !                  returns it as Y                               c
    !----------------------------------------------------------------c


    implicit none

    integer N,N_ARR,KLO,KHI,K,index
    real(8) XA(N_ARR),YA(N_ARR),Y2A(N_ARR),X,Y,H,A,B

    !       KLO=1
    !       KHI=N
    !  1    if (KHI-KLO.gt.1) then
    !          K=(KHI+KLO)/2
    !          if (XA(K).gt.X) then
    !             KHI=K
    !          else
    !             KLO=K
    !          endif
    !          goto 1
    !       endif
    !      print *,'inside subr: KLO=',KLO,', KHI=',KHI
    !
    KLO=index
    KHI=index+1
    !      print *,'before subr: KLO=',KLO,', KHI=',KHI
    !c      stop

    H=XA(KHI)-XA(KLO)
    A=(XA(KHI)-X)/H
    B=(X-XA(KLO))/H
    Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2) &
        /6.d0

end subroutine splint
subroutine svbksb(u,w,v,m,n,b, &
        tmp, &
        x)
    ! based on numerical recipes
    ! Marcel Sluiter, May 10 2002
    implicit none
    integer m,n,j
    real(8) b(m),u(m,n),v(n,n),w(n),x(n),tmp(n)
    intent(in) u,w,v,m,n,b
    intent(out) x,tmp

    tmp=0d0
    do j=1,n
        if(w(j).ne.0d0) tmp(j)=sum(u(:,j)*b(:))/w(j)
    enddo
    do j=1,n
        x(j)=sum(v(j,:)*tmp(:))
    enddo
end subroutine svbksb
subroutine svdcmp(m,n, &
        a, &
        tmp, &
        w,v)
    ! based on numerical recipes
    ! Marcel Sluiter, May 10 2002
    !      use m_inout
    implicit none
    integer m,n,i,its,j,jj,k,kk,nm
    real(8) a(m,n),v(n,n),w(n),anorm,c,f,g,h,s,scale,x,y,z,tmp(n), &
        pythag
    intent(in) m,n
    intent(out) w,v,tmp
    !i      include 'incmpi.h'

    g=0d0
    scale=0d0
    anorm=0d0
    do i=1,n
        kk=i+1
        tmp(i)=scale*g
        g=0d0
        s=0d0
        scale=0d0
        if(i.le.m)then
            do k=i,m
                scale=scale+abs(a(k,i))
            enddo
            if(scale.ne.0d0)then
                do k=i,m
                    a(k,i)=a(k,i)/scale
                    s=s+a(k,i)*a(k,i)
                enddo
                f=a(i,i)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,i)=f-g
                do j=kk,n
                    s=0d0
                    do k=i,m
                        s=s+a(k,i)*a(k,j)
                    enddo
                    f=s/h
                    do k=i,m
                        a(k,j)=a(k,j)+f*a(k,i)
                    enddo
                enddo
                do k=i,m
                    a(k,i)=scale*a(k,i)
                enddo
            endif
        endif
        w(i)=scale *g
        g=0d0
        s=0d0
        scale=0d0
        if((i.le.m).and.(i.ne.n))then
            do k=kk,n
                scale=scale+abs(a(i,k))
            enddo
            if(scale.ne.0d0)then
                do k=kk,n
                    a(i,k)=a(i,k)/scale
                    s=s+a(i,k)*a(i,k)
                enddo
                f=a(i,kk)
                g=-sign(sqrt(s),f)
                h=f*g-s
                a(i,kk)=f-g
                do k=kk,n
                    tmp(k)=a(i,k)/h
                enddo
                do j=kk,m
                    s=0d0
                    do k=kk,n
                        s=s+a(j,k)*a(i,k)
                    enddo
                    do k=kk,n
                        a(j,k)=a(j,k)+s*tmp(k)
                    enddo
                enddo
                do k=kk,n
                    a(i,k)=scale*a(i,k)
                enddo
            endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(tmp(i))))
    enddo
    do i=n,1,-1
        if(i.lt.n)then
            if(g.ne.0d0)then
                do j=kk,n
                    v(j,i)=(a(i,j)/a(i,kk))/g
                enddo
                do j=kk,n
                    s=0d0
                    do k=kk,n
                        s=s+a(i,k)*v(k,j)
                    enddo
                    do k=kk,n
                        v(k,j)=v(k,j)+s*v(k,i)
                    enddo
                enddo
            endif
            do j=kk,n
                v(i,j)=0d0
                v(j,i)=0d0
            enddo
        endif
        v(i,i)=1d0
        g=tmp(i)
        kk=i
    enddo
    do i=min(m,n),1,-1
        kk=i+1
        g=w(i)
        do j=kk,n
            a(i,j)=0d0
        enddo
        if(g.ne.0d0)then
            g=1d0/g
            do j=kk,n
                s=0d0
                do k=kk,m
                    s=s+a(k,i)*a(k,j)
                enddo
                f=(s/a(i,i))*g
                do k=i,m
                    a(k,j)=a(k,j)+f*a(k,i)
                enddo
            enddo
            do j=i,m
                a(j,i)=a(j,i)*g
            enddo
        else
            do j= i,m
                a(j,i)=0d0
            enddo
        endif
        a(i,i)=a(i,i)+1d0
    enddo
    do k=n,1,-1
        do its=1,30
            do kk=k,2,-1   !correction based on Kanno's
                nm=kk-1
                if((abs(tmp(kk))+anorm).eq.anorm)  goto 2
                if((abs(w(nm))+anorm).eq.anorm)  goto 1
            enddo
            kk=1   !correction based on Kanno's
            if((abs(tmp(kk))+anorm).eq.anorm)  goto 2   !correction
            1         c=0d0
            s=1d0
            do i=kk,k
                f=s*tmp(i)
                tmp(i)=c*tmp(i)
                if((abs(f)+anorm).eq.anorm) goto 2
                g=w(i)
                h=pythag(f,g)
                w(i)=h
                h=1d0/h
                c= (g*h)
                s=-(f*h)
                do j=1,m
                    y=a(j,nm)
                    z=a(j,i)
                    a(j,nm)=(y*c)+(z*s)
                    a(j,i)=-(y*s)+(z*c)
                enddo
            enddo
            2         z=w(k)
            if(kk.eq.k)then
                if(z.lt.0d0)then
                    w(k)=-z
                    v(:,k)=-v(:,k)
                endif
                goto 3
            endif
            !i          if(mype.eq.0.and.its.eq.30) write(error,'(a)')
            !i     &      '*WARNING* BROYDEN: no convergence in SVDCMP'
            x=w(kk)
            nm=k-1
            y=w(nm)
            g=tmp(nm)
            h=tmp(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2d0*h*y)
            g=pythag(f,1d0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=1d0
            s=1d0
            do j=kk,nm
                i=j+1
                g=tmp(i)
                y=w(i)
                h=s*g
                g=c*g
                z=pythag(f,h)
                tmp(j)=z
                c=f/z
                s=h/z
                f= (x*c)+(g*s)
                g=-(x*s)+(g*c)
                h=y*s
                y=y*c
                do jj=1,n
                    x=v(jj,j)
                    z=v(jj,i)
                    v(jj,j)= (x*c)+(z*s)
                    v(jj,i)=-(x*s)+(z*c)
                enddo
                z=pythag(f,h)
                w(j)=z
                if(z.ne.0d0)then
                    z=1d0/z
                    c=f*z
                    s=h*z
                endif
                f= (c*g)+(s*y)
                x=-(s*g)+(c*y)
                do jj=1,m
                    y=a(jj,j)
                    z=a(jj,i)
                    a(jj,j)= (y*c)+(z*s)
                    a(jj,i)=-(y*s)+(z*c)
                enddo
            enddo
            tmp(kk)=0d0
            tmp(k)=f
            w(k)=x
        enddo
        3       continue
    enddo
end subroutine svdcmp
subroutine svdinv(a,m,n,mdim,ndim, &
        ainv)
    ! computes the least squares inverse Ainv for a m*n matrix A
    ! method is based on Singular Value Decomposition as described in
    ! NUMERICAL RECIPES (2nd Ed) by W. Press et al. (Cambridge Univ. Press,
    ! New York, 1992), pp.51.
    ! the input is left intact
    ! A is defined as A(mdim,ndim), Ainv(ndim,mdim)
    ! Marcel Sluiter, Nov 28 2001

    implicit none
    real(8), parameter:: eps=1d-10  !cutoff for singularities
    integer mdim,ndim,i,n,m
    real(8) a(mdim,ndim),ainv(ndim,mdim)
    real(8), allocatable:: ac(:,:),y(:),wm(:),vm(:,:),tmp(:)
    intent(in) a,m,n,mdim,ndim
    intent(out) ainv

    allocate(ac(m,n),y(m),wm(n),vm(n,n),tmp(n))
    ac(1:m,1:n)=a(1:m,1:n)
    wm = 0d0
    tmp= 0d0
    vm = 0d0
    call svdcmp(m,n, &
        ac, &
        tmp, &
        wm,vm)
    where(wm(1:n).lt.eps) wm=0d0
    do i=1,m
        y=0d0
        y(i)=1d0
        call svbksb(ac,wm,vm,m,n,y, &
            tmp,ainv(i,1:n))
    enddo
    deallocate(ac,y,wm,vm,tmp)
end subroutine svdinv
      SUBROUTINE symnum(strng,type) 
      ! Subroutine to identify the type of "strng", numeric or symbolic 
      ! type returned is "INTG" for a string that contains a integer type 
      ! number (no decimal and no exponent marker); "NUME" if "strng"
      ! contains 
      ! a valid number real number; or "SYMB" otherwise. 
      ! valid characters in a number are: numbers from 0 to 9 
      ! the "+" and "-" symbols 
      ! decimal marker "." 
      ! "E" and "e" (exponent markers) 
      ! this routine expects input strings containing no leading blanks 
      ! Note: the presence of a double quote ('"') at the beginning of the 
      ! string is expected to be used to signal a string that would 
      ! otherwise be identified as being a scalar 
      !
      ! Vincent G (provided free on Yahoo.answers)
      
      implicit none

      CHARACTER*80 strng 
      CHARACTER*4 type 
      INTEGER lste,iper,i,jblnk
      
      iper=0 
      lste=0 
      jblnk=INDEX(strng,' ')-1 
      IF(jblnk.LT.0) jblnk=80 
      ! get first blank position: this is the depth of the input string 
      
      type='INTG' 
      ! assume this is an integer number 
      DO 10 i=1,jblnk 
      ! plain number detected 
      IF(strng(i:i).GE.'0' .AND. strng(i:i).LE.'9') GOTO 10 
      ! sign detection, must be before any number or right after an exponent
      ! marker 
      IF(strng(i:i).EQ.'+'.OR.strng(i:i).EQ.'-') THEN 
      ! sign after the exponent mark or sign in the first position 
      IF(lste.EQ.(i-1)) GOTO 10 
      ! else, not a number 
      GOTO 20 
      ! the letter "E" is present, may be text or maybe the exponent marker 
      ELSEIF(strng(i:i).EQ.'E'.OR.strng(i:i).EQ.'e') THEN 
      type='NUME' 
      ! downgrade to real type number 
      IF(lste.EQ.0) THEN 
      lste=i 
      GOTO 10 
      ELSE 
      ! there was another "E" before 
      GOTO 20 
      ENDIF 
      ! the period symbol is valid, but only once 
      ELSEIF(strng(i:i).EQ.'.') THEN 
      type='NUME' 
      ! downgrade to real type number 
      IF(iper.NE.0) GOTO 20 
      iper=i 
      ELSE 
      ! any other character is invalid for a number and is considered a
      ! separator 
      GOTO 20 
      ENDIF 
      10 CONTINUE 
      RETURN 
      20 type='SYMB' 
      RETURN 
      END 
subroutine testargs

    !-------------------------------------------------------------------c
    !
    !     Check if arguments have been supplied from command line.
    !
    !     Called by:     program MEAMfit
    !     Calls:         -
    !     Returns:       contjob, noOpt, readpotfile
    !
    !     Andrew Duff 2014
    !
    !-------------------------------------------------------------------c

    use m_generalinfo
    use m_optimization
    use m_meamparameters
    use m_filenames

    implicit none

    integer i,narg,num
    character*4 type
    character*20 name,name2
    character*80 filename

    !defaults:
    readpotfile=.false.
    contjob=.false.

    narg=command_argument_count()

    if (narg.gt.0) then
       i=1
       do
          call get_command_argument(i,name)
          select case(adjustl(name))
          case("-help","-h")
             i=i+1
             call get_command_argument(i,name2)
             if (adjustl(name2).eq."command-line") then
                print *
                print *,'Command-line options:'
                print *,'  -potfilein/-i filename   Read in potential from file: filename. Used either'
                print *,'                       with -onehot/-o (see below) to check fitted data, or'
                print *,'                       to provide starting point for further optimization.'
                print *,'  -noopt/-no           One iteration; no random para-gen (use in combination'
                print *,'                       with -potfilein/-i filename to check fitted data for a'
                print *,'                       given potential.'
             elseif (adjustl(name2).eq."settings-options") then
                print *,'Settings-file commands:'
                print *,'  Necessary to set:'
                print *,'  OPTPARAREADIN=1/2    The input method for specifying which parameters are'
                print *,'                       to be optimized. 1: select the parameters one-by-one'
                print *,'                       using command PARASTOOPT. 2: use keywords'
                print *,'                       (NTERMS_DENS, NTERMS_PP, NTERMS_EMB, etc) to specify'
                print *,'                       parameters.'
                print *,'  Optional:'
                print *,'  POTFILEIN=filename   Read in potential from file: filename.'
                print *,'  STARTPARAS=1/2       1: no random generation of potential parameters'
                print *,'                       (overrides all other commands relating to random'
                print *,'                       generation of potential parameters; this option requires'
                print *,'                       POTFILEIN). 2: enable random generation of potential'
                print *,'                       parameters (with the particular parameters to optimized'
                print *,'                       determined by other commands in settings files)'
                print *,'  NOOPT=TRUE           No optimization'
             else
                print *
                print *,'Options:'
                print *,'  -help command-line   Show possible command-line options.'
                print *,'  -help settings-options  Show possible variables for inclusion in settings'
                print *,'                       file.'
             endif
             stop
          case("-noopt","-no")
             noOpt=.true.
         !case("-ncutoff","-nc")
         !   i=i+1
         !   call get_command_argument(i,name)
         !   call symnum(name,type)
         !   if (type.ne."INTG") then
         !      print *,'Warning, no integer after -ncutoff/-nc, STOPPING.'
         !      stop
         !   endif
         !   read(name, '(i10)' ) nCutoffOverride
          case("-potfilein","-i")
             i=i+1
             call get_command_argument(i,startparameterfile)
             call symnum(startparameterfile,type)
             if (type.ne."SYMB") then
                print *,'Warning, expected filename after -potin/-i, STOPPING.'
                stop
             endif
             readpotfile=.true.
             !print *,'file=',startparameterfile
             !stop
          case("-cont")
             contjob=.true.
          case default
             print *,'Command-line option un-known; stopping.'
             stop
          end select
          if (i.eq.narg) then
             exit
          endif
          i=i+1
       enddo
    endif

end subroutine testargs
      subroutine assst ( iv, liv, lv, v )
      
      !*****************************************************************************80
      !
      !! ASSST assesses a candidate step.
      !
      !  Discussion:
      !
      !    This subroutine is called by an unconstrained minimization
      !    routine to assess the next candidate step.  it may recommend one
      !    of several courses of action, such as accepting the step, recom-
      !    puting it using the same or a new quadratic model, or halting due
      !    to convergence or false convergence.  See the return code listing
      !    below.
      !
      !  Reference:
      !
      !    John Dennis, David Gay, Roy Welsch,
      !    An Adaptive Nonlinear Least-squares Algorithm,
      !    ACM Transactions on Mathematical Software,
      !    Volume 7, Number 3, 1981.
      !
      !    M J D Powell,
      !    A Fortran Subroutine for Solving Systems of Nonlinear Algebraic
      !    Equations,
      !    in Numerical Methods for Nonlinear Algebraic Equations,
      !    edited by Philip Rabinowitz,
      !    Gordon and Breach, London, 1970.
      !
      !  Parameters:
      !
      !  iv (i/o) integer parameter and scratch vector -- see description
      !             below of iv values referenced.
      !
      ! liv (in)  length of iv array.
      !
      !  lv (in)  length of v array.
      !
      !   v (i/o) real parameter and scratch vector -- see description
      !             below of v values referenced.
      !
      !   iv values referenced
      !
      !    iv(irc) (i/o) on input for the first step tried in a new iteration,
      !             iv(irc) should be set to 3 or 4 (the value to which it is
      !             set when step is definitely to be accepted).  on input
      !             after step has been recomputed, iv(irc) should be
      !             unchanged since the previous return of assst.
      !                on output, iv(irc) is a return code having one of the
      !             following values...
      !                  1 = switch models or try smaller step.
      !                  2 = switch models or accept step.
      !                  3 = accept step and determine v(radfac) by gradient
      !                       tests.
      !                  4 = accept step, v(radfac) has been determined.
      !                  5 = recompute step (using the same model).
      !                  6 = recompute step with radius = v(lmaxs) but do not
      !                       evaulate the objective function.
      !                  7 = x-convergence (see v(xctol)).
      !                  8 = relative function convergence (see v(rfctol)).
      !                  9 = both x- and relative function convergence.
      !                 10 = absolute function convergence (see v(afctol)).
      !                 11 = singular convergence (see v(lmaxs)).
      !                 12 = false convergence (see v(xftol)).
      !                 13 = iv(irc) was out of range on input.
      !             return code i has precdence over i+1 for i = 9, 10, 11.
      ! iv(mlstgd) (i/o) saved value of iv(model).
      !  iv(model) (i/o) on input, iv(model) should be an integer identifying
      !             the current quadratic model of the objective function.
      !             if a previous step yielded a better function reduction,
      !             then iv(model) will be set to iv(mlstgd) on output.
      ! iv(nfcall) (in)  invocation count for the objective function.
      ! iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
      !             function reduction this iteration.  iv(nfgcal) remains
      !             unchanged until a function reduction is obtained.
      ! iv(radinc) (i/o) the number of radius increases (or minus the number
      !             of decreases) so far this iteration.
      ! iv(restor) (out) set to 1 if v(f) has been restored and x should be
      !             restored to its initial value, to 2 if x should be saved,
      !             to 3 if x should be restored from the saved value, and to
      !             0 otherwise.
      !  iv(stage) (i/o) count of the number of models tried so far in the
      !             current iteration.
      ! iv(stglim) (in)  maximum number of models to consider.
      ! iv(switch) (out) set to 0 unless a new model is being tried and it
      !             gives a smaller function value than the previous model,
      !             in which case assst sets iv(switch) = 1.
      ! iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
      !             overflow).
      !   iv(xirc) (i/o) value that iv(irc) would have in the absence of
      !             convergence, false convergence, and oversized steps.
      !
      !   v values referenced
      !
      ! v(afctol) (in)  absolute function convergence tolerance.  if the
      !             absolute value of the current function value v(f) is less
      !             than v(afctol), then assst returns with iv(irc) = 10.
      ! v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
      !             nonzero.
      ! v(dstnrm) (in)  the 2-norm of d*step.
      ! v(dstsav) (i/o) value of v(dstnrm) on saved step.
      !   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
      !             i.e., for v(nreduc) >= 0).
      !      v(f) (i/o) on both input and output, v(f) is the objective func-
      !             tion value at x.  if x is restored to a previous value,
      !             then v(f) is restored to the corresponding value.
      !   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
      !             value of v(f) if an earlier step gave a bigger function
      !             decrease, and for the input value of v(f) otherwise).
      ! v(flstgd) (i/o) saved value of v(f).
      !     v(f0) (in)  objective function value at start of iteration.
      ! v(gtslst) (i/o) value of v(gtstep) on saved step.
      ! v(gtstep) (in)  inner product between step and gradient.
      ! v(incfac) (in)  minimum factor by which to increase radius.
      !  v(lmaxs) (in)  maximum reasonable step size (and initial step bound).
      !             if the actual function decrease is no more than twice
      !             what was predicted, if a return with iv(irc) = 7, 8, 9,
      !             or 10 does not occur, if v(dstnrm) > v(lmaxs), and if
      !             v(preduc) <= v(sctol) * abs(v(f0)), then assst re-
      !             turns with iv(irc) = 11.  if so doing appears worthwhile,
      !             then assst repeats this test with v(preduc) computed for
      !             a step of length v(lmaxs) (by a return with iv(irc) = 6).
      ! v(nreduc) (i/o)  function reduction predicted by quadratic model for
      !             newton step.  if assst is called with iv(irc) = 6, i.e.,
      !             if v(preduc) has been computed with radius = v(lmaxs) for
      !             use in the singular convervence test, then v(nreduc) is
      !             set to -v(preduc) before the latter is restored.
      ! v(plstgd) (i/o) value of v(preduc) on saved step.
      ! v(preduc) (i/o) function reduction predicted by quadratic model for
      !             current step.
      ! v(radfac) (out) factor to be used in determining the new radius,
      !             which should be v(radfac)*dst, where  dst  is either the
      !             output value of v(dstnrm) or the 2-norm of
      !             diag(newd)*step  for the output value of step and the
      !             updated version, newd, of the scale vector d.  for
      !             iv(irc) = 3, v(radfac) = 1.0D+00 is returned.
      ! v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
      !             value of v(dstnrm) -- suggested value = 0.1.
      ! v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
      !  v(reldx) (in) scaled relative change in x caused by step, computed
      !             (e.g.) by function  reldst  as
      !                 max (d(i)*abs(x(i)-x0(i)), 1 <= i <= p) /
      !                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p).
      ! v(rfctol) (in)  relative function convergence tolerance.  if the
      !             actual function reduction is at most twice what was pre-
      !             dicted and  v(nreduc) <= v(rfctol)*abs(v(f0)),  then
      !             assst returns with iv(irc) = 8 or 9.
      ! v(stppar) (in)  marquardt parameter -- 0 means full newton step.
      ! v(tuner1) (in)  tuning constant used to decide if the function
      !             reduction was much less than expected.  suggested
      !             value = 0.1.
      ! v(tuner2) (in)  tuning constant used to decide if the function
      !             reduction was large enough to accept step.  suggested
      !             value = 10**-4.
      ! v(tuner3) (in)  tuning constant used to decide if the radius
      !             should be increased.  suggested value = 0.75.
      !  v(xctol) (in)  x-convergence criterion.  if step is a newton step
      !             (v(stppar) = 0) having v(reldx) <= v(xctol) and giving
      !             at most twice the predicted function decrease, then
      !             assst returns iv(irc) = 7 or 9.
      !  v(xftol) (in)  false convergence tolerance.  if step gave no or only
      !             a small function decrease and v(reldx) <= v(xftol),
      !             then assst returns with iv(irc) = 12.
      !
      !  notes
      !
      !   application and usage restrictions
      !
      !        this routine is called as part of the nl2sol (nonlinear
      !     least-squares) package.  it may be used in any unconstrained
      !     minimization solver that uses dogleg, goldfeld-quandt-trotter,
      !     or levenberg-marquardt steps.
      !
      !   algorithm notes
      !
      !        see (1) for further discussion of the assessing and model
      !     switching strategies.  while nl2sol considers only two models,
      !     assst is designed to handle any number of models.
      !
      !   usage notes
      !
      !        on the first call of an iteration, only the i/o variables
      !     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
      !     v(preduc) need have been initialized.  between calls, no i/o
      !     values execpt step, x, iv(model), v(f) and the stopping toler-
      !     ances should be changed.
      !        after a return for convergence or false convergence, one can
      !     change the stopping tolerances and call assst again, in which
      !     case the stopping tests will be repeated.
      !
      !   history
      !
      !        john dennis designed much of this routine, starting with
      !     ideas in (2). roy welsch suggested the model switching strategy.
      !        david gay and stephen peters cast this subroutine into a more
      !     portable form (winter 1977), and david gay cast it into its
      !     present form (fall 1978).
      !
        integer liv
        integer lv
      
        integer iv(liv)
        real ( kind = 8 ) v(lv)
        logical goodx
        integer i, nfc
        real ( kind = 8 ) emax, emaxs, gts, rfac1, xmax
        real ( kind = 8 ) half, one, onep2, two
        integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0
        integer gtslst, gtstep, incfac, irc, lmaxs, mlstgd, model, nfcall
        integer nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn
        integer rdfcmx, reldx, restor, rfctol, sctol, stage, stglim
        integer stppar, switch, toobig, tuner1, tuner2, tuner3, xctol
        integer xftol, xirc
      
        parameter ( half=0.5d+0, one=1.d+0, onep2=1.2d+0, two=2.d+0)
        parameter ( irc=29, mlstgd=32, model=5, nfcall=6, nfgcal=7 )
        parameter ( radinc=8, restor=9, stage=10, stglim=11, switch=12 )
        parameter ( toobig=2, xirc=13)
        parameter (afctol=31, decfac=22, dstnrm=2, dst0=3, dstsav=18 )
        parameter (f=10, fdif=11, flstgd=12, f0=13, gtslst=14, gtstep=4 )
        parameter (incfac=23, lmaxs=36, nreduc=6, plstgd=15, preduc=7 )
        parameter (radfac=16, rdfcmn=24, rdfcmx=25, reldx=17, rfctol=32 )
        parameter (sctol=37, stppar=5, tuner1=26, tuner2=27, tuner3=28 )
        parameter (xctol=33, xftol=34)
      
        nfc = iv(nfcall)
        iv(switch) = 0
        iv(restor) = 0
        rfac1 = one
        goodx = .true.
        i = iv(irc)
      
        if (i >= 1 .and. i <= 12) then
              go to (20,30,10,10,40,280,220,220,220,220,220,170), i
        end if
      
        iv(irc) = 13
        return
      !
      !  Initialize for new iteration.
      !
       10   iv(stage) = 1
        iv(radinc) = 0
        v(flstgd) = v(f0)
        if (iv(toobig) == 0) go to 110
           iv(stage) = -1
           iv(xirc) = i
           go to 60
      !
      !  Step was recomputed with new model or smaller radius
      !  first decide which
      !
       20   if (iv(model) /= iv(mlstgd)) go to 30
      !
      !  Old model retained, smaller radius tried
      !  do not consider any more new models this iteration
      !
           iv(stage) = iv(stglim)
           iv(radinc) = -1
           go to 110
      !
      !  A new model is being tried.  decide whether to keep it.
      !
       30   iv(stage) = iv(stage) + 1
      !
      !  Now we add the possibiltiy that step was recomputed with
      !  the same model, perhaps because of an oversized step.
      !
       40   if (iv(stage) > 0) go to 50
      !
      !  Step was recomputed because it was too big.
      !
           if (iv(toobig) /= 0) go to 60
      !
      !  Restore iv(stage) and pick up where we left off.
      !
           iv(stage) = -iv(stage)
           i = iv(xirc)
           go to (20, 30, 110, 110, 70), i
      
       50   if (iv(toobig) == 0) go to 70
      !
      !  Handle oversize step
      !
        if (iv(radinc) > 0) go to 80
           iv(stage) = -iv(stage)
           iv(xirc) = iv(irc)
      
       60      v(radfac) = v(decfac)
           iv(radinc) = iv(radinc) - 1
           iv(irc) = 5
           iv(restor) = 1
           return
      
       70   if (v(f) < v(flstgd)) go to 110
      !
      !  The new step is a loser.  restore old model.
      !
        if (iv(model) == iv(mlstgd)) go to 80
           iv(model) = iv(mlstgd)
           iv(switch) = 1
      !
      !  Restore step, etc. only if a previous step decreased v(f).
      !
       80   if (v(flstgd) >= v(f0)) go to 110
           iv(restor) = 1
           v(f) = v(flstgd)
           v(preduc) = v(plstgd)
           v(gtstep) = v(gtslst)
           if (iv(switch) == 0) rfac1 = v(dstnrm) / v(dstsav)
           v(dstnrm) = v(dstsav)
           nfc = iv(nfgcal)
           goodx = .false.
      
       110  v(fdif) = v(f0) - v(f)
        if (v(fdif) > v(tuner2) * v(preduc)) go to 140
        if(iv(radinc)>0) go to 140
      !
      !         no (or only a trivial) function decrease
      !         so try new model or smaller radius
      !
           if (v(f) < v(f0)) go to 120
                iv(mlstgd) = iv(model)
                v(flstgd) = v(f)
                v(f) = v(f0)
                iv(restor) = 1
                go to 130
       120     iv(nfgcal) = nfc
       130     iv(irc) = 1
           if (iv(stage) < iv(stglim)) go to 160
                iv(irc) = 5
                iv(radinc) = iv(radinc) - 1
                go to 160
      !
      !  Nontrivial function decrease achieved
      !
       140  iv(nfgcal) = nfc
        rfac1 = 1.0D+00
        v(dstsav) = v(dstnrm)
        if (v(fdif) > v(preduc)*v(tuner1)) go to 190
      !
      !  Decrease was much less than predicted -- either change models
      !  or accept step with decreased radius.
      !
        if (iv(stage) >= iv(stglim)) go to 150
      !
      !  Consider switching models
      !
           iv(irc) = 2
           go to 160
      !
      !  Accept step with decreased radius
      !
       150  iv(irc) = 4
      !
      !   set v(radfac) to fletcher*s decrease factor
      !
       160  iv(xirc) = iv(irc)
        emax = v(gtstep) + v(fdif)
        v(radfac) = half * rfac1
      
        if (emax < v(gtstep)) then
          v(radfac) = rfac1 * max (v(rdfcmn),half * v(gtstep)/emax)
        end if
      !
      !  Do false convergence test
      !
       170  if (v(reldx) <= v(xftol)) go to 180
           iv(irc) = iv(xirc)
           if (v(f) < v(f0)) go to 200
                go to 230
      
       180  iv(irc) = 12
        go to 240
      !
      !  Handle good function decrease
      !
       190  if (v(fdif) < (-v(tuner3) * v(gtstep))) go to 210
      !
      !  Increasing radius looks worthwhile.  see if we just
      !  recomputed step with a decreased radius or restored step
      !  after recomputing it with a larger radius.
      !
        if (iv(radinc) < 0) go to 210
        if (iv(restor) == 1) go to 210
      !
      !  We did not.  try a longer step unless this was a newton step.
      !
           v(radfac) = v(rdfcmx)
           gts = v(gtstep)
           if (v(fdif) < (half/v(radfac) - 1.0D+00 ) * gts) then
             v(radfac) = max (v(incfac), half*gts/(gts + v(fdif)))
           end if
           iv(irc) = 4
           if (v(stppar) == 0.0D+00 ) go to 230
           if (v(dst0) >= 0.0D+00 .and. (v(dst0) < two*v(dstnrm) &
                    .or. v(nreduc) < onep2*v(fdif)))  then
             go to 230
           end if
      !
      !  Step was not a newton step.  recompute it with a larger radius.
      !
                iv(irc) = 5
                iv(radinc) = iv(radinc) + 1
      !
      !  Save values corresponding to good step
      !
       200  v(flstgd) = v(f)
        iv(mlstgd) = iv(model)
        if (iv(restor) /= 1) iv(restor) = 2
        v(dstsav) = v(dstnrm)
        iv(nfgcal) = nfc
        v(plstgd) = v(preduc)
        v(gtslst) = v(gtstep)
        go to 230
      !
      !  Accept step with radius unchanged.
      !
       210  v(radfac) = 1.0D+00
        iv(irc) = 3
        go to 230
      !
      !  Come here for a restart after convergence.
      !
       220  iv(irc) = iv(xirc)
        if (v(dstsav) >= 0.0D+00 ) go to 240
           iv(irc) = 12
           go to 240
      !
      !  Perform convergence tests.
      !
       230  iv(xirc) = iv(irc)
       240  if (iv(restor) == 1 .and. v(flstgd) < v(f0)) iv(restor) = 3
        if (abs(v(f)) < v(afctol)) iv(irc) = 10
      
        if (half * v(fdif) > v(preduc)) then
          return
        end if
      
        emax = v(rfctol) * abs(v(f0))
        emaxs = v(sctol) * abs(v(f0))
        if (v(dstnrm) > v(lmaxs) .and. v(preduc) <= emaxs) then
          iv(irc) = 11
        end if
        if (v(dst0) < 0.0D+00 ) go to 250
        i = 0
      
        if ((v(nreduc) > 0.0D+00 .and. v(nreduc) <= emax) .or. &
            (v(nreduc) == 0.0D+00 .and. v(preduc) == 0.0D+00 )) then
          i = 2
        end if
      
        if (v(stppar) == 0.0D+00 .and. v(reldx) <= v(xctol) .and. goodx) then
          i = i + 1
        end if
      
        if (i > 0) iv(irc) = i + 6
      !
      !  Consider recomputing step of length v(lmaxs) for singular
      !  convergence test.
      !
       250  if (iv(irc) > 5 .and. iv(irc) /= 12) then
           return
        end if
      
        if (v(dstnrm) > v(lmaxs)) go to 260
           if (v(preduc) >= emaxs) then
             return
           end if
                if (v(dst0) <= 0.0D+00 ) go to 270
                     if (half * v(dst0) <= v(lmaxs)) then
                       return
                     end if
                          go to 270
       260  if (half * v(dstnrm) <= v(lmaxs)) then
              return
            end if
        xmax = v(lmaxs) / v(dstnrm)
        if (xmax * (two - xmax) * v(preduc) >= emaxs) then
          return
        end if
       270  if (v(nreduc) < 0.0D+00 ) go to 290
      !
      !   recompute v(preduc) for use in singular convergence test
      !
        v(gtslst) = v(gtstep)
        v(dstsav) = v(dstnrm)
        if (iv(irc) == 12) v(dstsav) = -v(dstsav)
        v(plstgd) = v(preduc)
        i = iv(restor)
        iv(restor) = 2
        if (i == 3) iv(restor) = 0
        iv(irc) = 6
        return
      !
      !  Perform singular convergence test with recomputed v(preduc)
      !
       280  v(gtstep) = v(gtslst)
        v(dstnrm) = abs(v(dstsav))
        iv(irc) = iv(xirc)
        if (v(dstsav) <= 0.0D+00 ) iv(irc) = 12
        v(nreduc) = -v(preduc)
        v(preduc) = v(plstgd)
        iv(restor) = 3
      
       290  if (-v(nreduc) <= v(rfctol) * abs(v(f0))) iv(irc) = 11
      
        return
      end
      subroutine dbdog ( dig, lv, n, nwtstp, step, v )
      
      !*****************************************************************************80
      !
      !! DBDOG: compute a double dogleg step.
      !
      !  Discussion:
      !
      !    This subroutine computes a candidate step (for use in an
      !    unconstrained minimization code) by the double dogleg algorithm of
      !    dennis and mei (ref. 1), which is a variation on powell*s dogleg
      !    scheme (ref. 2, p. 95).
      !
      !    let  g  and  h  be the current gradient and hessian approxima-
      !    tion respectively and let d be the current scale vector.  this
      !    routine assumes dig = diag(d)**-2 * g  and  nwtstp = h**-1 * g.
      !    the step computed is the same one would get by replacing g and h
      !    by  diag(d)**-1 * g  and  diag(d)**-1 * h * diag(d)**-1,
      !    computing step, and translating step back to the original
      !    variables, i.e., premultiplying it by diag(d)**-1.
      !
      !  Reference:
      !
      !    John Dennis, Howell Mei,
      !    Two New Unconstrained Optimization Algorithms Which Use
      !    Function and Gradient Values,
      !    Journal of Optimization Theory and Applications,
      !    Volume 28, pages 453-482, 1979.
      !
      !    M J D Powell,
      !    A Hybrid Method for Non-linear Equations,
      !    in Numerical Methods for Non-linear Equations,
      !    edited by Philip Rabinowitz,
      !    Gordon and Breach, London, 1970.
      !
      !  Parameters:
      !
      !    dig (input) diag(d)**-2 * g -- see algorithm notes.
      !      g (input) the current gradient vector.
      !     lv (input) length of v.
      !      n (input) number of components in  dig, g, nwtstp,  and  step.
      ! nwtstp (input) negative newton step -- see algorithm notes.
      !   step (output) the computed step.
      !      v (i/o) values array, the following components of which are
      !             used here...
      ! v(bias)   (input) bias for relaxed newton step, which is v(bias) of
      !             the way from the full newton to the fully relaxed newton
      !             step.  recommended value = 0.8 .
      ! v(dgnorm) (input) 2-norm of diag(d)**-1 * g -- see algorithm notes.
      ! v(dstnrm) (output) 2-norm of diag(d) * step, which is v(radius)
      !             unless v(stppar) = 0 -- see algorithm notes.
      ! v(dst0) (input) 2-norm of diag(d) * nwtstp -- see algorithm notes.
      ! v(grdfac) (output) the coefficient of  dig  in the step returned --
      !             step(i) = v(grdfac)*dig(i) + v(nwtfac)*nwtstp(i).
      ! v(gthg)   (input) square-root of (dig**t) * (hessian) * dig -- see
      !             algorithm notes.
      ! v(gtstep) (output) inner product between g and step.
      ! v(nreduc) (output) function reduction predicted for the full newton
      !             step.
      ! v(nwtfac) (output) the coefficient of  nwtstp  in the step returned --
      !             see v(grdfac) above.
      ! v(preduc) (output) function reduction predicted for the step returned.
      ! v(radius) (input) the trust region radius.  d times the step returned
      !             has 2-norm v(radius) unless v(stppar) = 0.
      ! v(stppar) (output) code telling how step was computed... 0 means a
      !             full newton step.  between 0 and 1 means v(stppar) of the
      !             way from the newton to the relaxed newton step.  between
      !             1 and 2 means a true double dogleg step, v(stppar) - 1 of
      !             the way from the relaxed newton to the Cauchy step.
      !             greater than 2 means 1 / (v(stppar) - 1) times the Cauchy
      !             step.
      !
        integer lv
        integer n
      
        real ( kind = 8 ) dig(n), nwtstp(n), step(n), v(lv)
        external dotprd, v2norm
        real ( kind = 8 ) dotprd, v2norm
        real ( kind = 8 ) cfact, cnorm, ctrnwt, ghinvg, femnsq, gnorm
        real ( kind = 8 ) nwtnrm, relax, rlambd, t, t1, t2
        real ( kind = 8 ) half, two
        integer bias, dgnorm, dstnrm, dst0, grdfac, gthg, gtstep
        integer nreduc, nwtfac, preduc, radius, stppar
        parameter (half=0.5d+0, two=2.d+0)
        parameter (bias=43, dgnorm=1, dstnrm=2, dst0=3, grdfac=45 )
        parameter ( gthg=44, gtstep=4, nreduc=6, nwtfac=46, preduc=7 )
        parameter ( radius=8, stppar=5)
      
        nwtnrm = v(dst0)
        rlambd = 1.0D+00
        if (nwtnrm > 0.0D+00 ) rlambd = v(radius) / nwtnrm
        gnorm = v(dgnorm)
        ghinvg = two * v(nreduc)
        v(grdfac) = 0.0D+00
        v(nwtfac) = 0.0D+00
        if (rlambd < 1.0D+00 ) go to 30
      !
      !  The Newton step is inside the trust region.
      !
           v(stppar) = 0.0D+00
           v(dstnrm) = nwtnrm
           v(gtstep) = -ghinvg
           v(preduc) = v(nreduc)
           v(nwtfac) = -1.0D+00
           step(1:n) = -nwtstp(1:n)
           return
      
       30   v(dstnrm) = v(radius)
        cfact = (gnorm / v(gthg))**2
      !
      !  Cauchy step = -cfact * g.
      !
        cnorm = gnorm * cfact
        relax = 1.0D+00 - v(bias) * ( 1.0D+00 - gnorm*cnorm/ghinvg)
        if (rlambd < relax) go to 50
      !
      !  Step is between relaxed Newton and full Newton steps.
      !
           v(stppar) = 1.0D+00 -  (rlambd - relax) / ( 1.0D+00 - relax)
           t = -rlambd
           v(gtstep) = t * ghinvg
           v(preduc) = rlambd * ( 1.0D+00 - half*rlambd) * ghinvg
           v(nwtfac) = t
           step(1:n) = t * nwtstp(1:n)
           return
      
       50   if (cnorm < v(radius)) go to 70
      !
      !  The Cauchy step lies outside the trust region --
      !  step = scaled Cauchy step.
      !
           t = -v(radius) / gnorm
           v(grdfac) = t
           v(stppar) = 1.0D+00  +  cnorm / v(radius)
           v(gtstep) = -v(radius) * gnorm
        v(preduc) = v(radius)*(gnorm - half*v(radius)*(v(gthg)/gnorm)**2)
           step(1:n) = t * dig(1:n)
           return
      !
      !  Compute dogleg step between Cauchy and relaxed Newton
      !  femur = relaxed newton step minus Cauchy step.
      !
       70   ctrnwt = cfact * relax * ghinvg / gnorm
      !
      !  ctrnwt = inner product of Cauchy and relaxed Newton steps,
      !  scaled by gnorm**-1.
      !
        t1 = ctrnwt - gnorm*cfact**2
      !
      !  t1 = inner prod. of femur and Cauchy step, scaled by gnorm**-1.
      !
        t2 = v(radius)*(v(radius)/gnorm) - gnorm*cfact**2
        t = relax * nwtnrm
        femnsq = (t/gnorm)*t - ctrnwt - t1
      !
      !  femnsq = square of 2-norm of femur, scaled by gnorm**-1.
      !
        t = t2 / (t1 + sqrt(t1**2 + femnsq*t2))
      !
      !  Dogleg step  =  Cauchy step  +  t * femur.
      !
        t1 = (t - 1.0D+00 ) * cfact
        v(grdfac) = t1
        t2 = -t * relax
        v(nwtfac) = t2
        v(stppar) = two - t
        v(gtstep) = t1*gnorm**2 + t2*ghinvg
        v(preduc) = -t1*gnorm * ((t2 + 1.0D+00 )*gnorm) &
                        - t2 * ( 1.0D+00 + half*t2)*ghinvg &
                         - half * (v(gthg)*t1)**2
      
        step(1:n) = t1 * dig(1:n) + t2 * nwtstp(1:n)
      
        return
      end
      subroutine deflt ( alg, iv, liv, lv, v )
      
      !*****************************************************************************80
      !
      !! DEFLT: supply default values to IV and V.
      !
      !  Discussion:
      !
      !   ALG = 1 means regression constants.
      !   ALG = 2 means general unconstrained optimization constants.
      !
        integer liv
        integer lv
      
        integer alg
        integer iv(liv)
        real ( kind = 8 ) v(lv)
        external vdflt
        integer miv, mv
        integer miniv(2), minv(2)
        integer algsav, covprt, covreq, dtype, hc, ierr, inith, inits
        integer ipivot, ivneed, lastiv, lastv, lmat, mxfcal, mxiter
        integer nfcov, ngcov, nvdflt, outlev, parprt, parsav, perm
        integer prunit, qrtyp, rdreq, rmat, solprt, statpr, vneed
        integer vsave, x0prt
      
        parameter (algsav=51, covprt=14, covreq=15, dtype=16, hc=71 )
        parameter (ierr=75, inith=25, inits=25, ipivot=76, ivneed=3 )
        parameter (lastiv=44, lastv=45, lmat=42, mxfcal=17, mxiter=18 )
        parameter (nfcov=52, ngcov=53, nvdflt=50, outlev=19, parprt=20 )
        parameter (parsav=49, perm=58, prunit=21, qrtyp=80, rdreq=57 )
        parameter (rmat=78, solprt=22, statpr=23, vneed=4, vsave=60 )
        parameter (x0prt=24)
      
        data miniv(1)/80/, miniv(2)/59/, minv(1)/98/, minv(2)/71/
      
        if ( alg < 1 .or. 2 < alg ) then
          iv(1) = 67
          return
        end if
      
        miv = miniv(alg)
      
        if ( liv < miv ) then
          iv(1) = 15
          return
        end if
      
        mv = minv(alg)
      
        if ( lv < mv ) then
          iv(1) = 16
          return
        end if
      
        call vdflt(alg, lv, v)
        iv(1) = 12
        iv(algsav) = alg
        iv(ivneed) = 0
        iv(lastiv) = miv
        iv(lastv) = mv
        iv(lmat) = mv + 1
        iv(mxfcal) = 200
        iv(mxiter) = 150
        iv(outlev) = 1
        iv(parprt) = 1
        iv(perm) = miv + 1
        iv(prunit) = 6
        iv(solprt) = 1
        iv(statpr) = 1
        iv(vneed) = 0
        iv(x0prt) = 1
      !
      !  General optimization values.
      !
        if ( 2 <= alg ) then
      
          iv(dtype) = 0
          iv(inith) = 1
          iv(nfcov) = 0
          iv(ngcov) = 0
          iv(nvdflt) = 25
          iv(parsav) = 47
      !
      !  Regression values.
      !
        else
      
          iv(covprt) = 3
          iv(covreq) = 1
          iv(dtype) = 1
          iv(hc) = 0
          iv(ierr) = 0
          iv(inits) = 0
          iv(ipivot) = 0
          iv(nvdflt) = 32
          iv(parsav) = 67
          iv(qrtyp) = 1
          iv(rdreq) = 3
          iv(rmat) = 0
          iv(vsave) = 58
      
        end if
      
        return
      end
      function dotprd ( p, x, y )
      
      !*****************************************************************************80
      !
      !! DOTPRD returns the inner product of vectors X and Y.
      !
        integer p
      
        real ( kind = 8 ) dotprd
        integer i
        real ( kind = 8 ) rmdcon
        real ( kind = 8 ), save :: sqteta = 0.0D+00
        real ( kind = 8 ) t
        real ( kind = 8 ) x(p)
        real ( kind = 8 ) y(p)
      
        dotprd = 0.0D+00
      
        if ( sqteta == 0.0D+00 ) then
          sqteta = rmdcon(2)
        end if
      
        do i = 1, p
      
          t = max ( abs ( x(i) ), abs ( y(i) ) )
          if ( t > 1.0D+00 ) go to 10
          if (t < sqteta) go to 20
          t = (x(i)/sqteta)*y(i)
          if (abs(t) < sqteta) go to 20
       10   dotprd = dotprd + x(i)*y(i)
      
       20 continue
      
        end do
      
        return
      end
      subroutine dupdu ( d, hdiag, iv, liv, lv, n, v )
      
      !*****************************************************************************80
      !
      !! DUPDU: update scale vector D for HUMSL.
      !
      !  Modified:
      !
      !    20 February 2006
      !
        integer liv
        integer lv
        integer n
      
        real ( kind = 8 ) d(n)
        integer d0i
        integer, parameter :: dfac = 41
        integer, parameter :: dtol = 59
        integer, parameter :: dtype = 16
        integer dtoli
        real ( kind = 8 ) hdiag(n)
        integer i
        integer iv(liv)
        integer, parameter :: niter = 31
        real ( kind = 8 ) t
        real ( kind = 8 ) v(lv)
        real ( kind = 8 ) vdfac
      
        i = iv(dtype)
      
        if ( i /= 1 ) then
          if ( 0 < iv(niter) ) then
            return
          end if
        end if
      
        dtoli = iv(dtol)
        d0i = dtoli + n
        vdfac = v(dfac)
      
        do i = 1, n
      
          t = max ( sqrt ( abs ( hdiag(i) ) ), vdfac * d(i) )
      
          if ( t < v(dtoli) ) then
            t = max ( v(dtoli), v(d0i) )
          end if
      
          d(i) = t
          dtoli = dtoli + 1
          d0i = d0i + 1
      
        end do
      
        return
      end
      subroutine gqtst ( d, dig, dihdi, ka, l, p, step, v, w )
      
      !*****************************************************************************80
      !
      !! GQTST: compute Goldfeld-Quandt-Trotter step by More-Hebden technique.
      !
      !  Discussion:
      !
      !    Given the (compactly stored) lower triangle of a scaled
      !    hessian (approximation) and a nonzero scaled gradient vector,
      !    this subroutine computes a goldfeld-quandt-trotter step of
      !    approximate length v(radius) by the more-hebden technique.  in
      !    other words, step is computed to (approximately) minimize
      !    psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the
      !    2-norm of d*step is at most (approximately) v(radius), where
      !    g  is the gradient,  h  is the hessian, and  d  is a diagonal
      !    scale matrix whose diagonal is stored in the parameter d.
      !    (gqtst assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.)
      !
      !    the desired g-q-t step (ref. 2, 3, 4, 6) satisfies
      !    (h + alpha*d**2)*step = -g  for some nonnegative alpha such that
      !    h + alpha*d**2 is positive semidefinite.  alpha and step are
      !    computed by a scheme analogous to the one described in ref. 5.
      !    estimates of the smallest and largest eigenvalues of the hessian
      !    are obtained from the gerschgorin circle theorem enhanced by a
      !    simple form of the scaling described in ref. 7.  cases in which
      !    h + alpha*d**2 is nearly (or exactly) singular are handled by
      !    the technique discussed in ref. 2.  in these cases, a step of
      !    (exact) length v(radius) is returned for which psi(step) exceeds
      !    its optimal value by less than -v(epslon)*psi(step).  the test
      !    suggested in ref. 6 for detecting the special case is performed
      !    once two matrix factorizations have been done -- doing so sooner
      !    seems to degrade the performance of optimization routines that
      !    call this routine.
      !
      !  Reference:
      !
      !    John Dennis, David Gay, Roy Welsch,
      !    An Adaptive Nonlinear Least-squares Algorithm,
      !    ACM Transactions on Mathematical Software,
      !    Volume 7, Number 3, 1981.
      !
      !    David Gay,
      !    Computing Optimal Locally Constrained Steps,
      !    SIAM Journal on Scientific and Statistical Computing,
      !    Volume 2, pages 186-197, 1981.
      !
      !    S M Goldfeld, R E Quandt, H F Trotter,
      !    Maximization by Quadratic Hill-climbing,
      !    Econometrica,
      !    Volume 34, pages 541-551, 1966.
      !
      !    M D Hebden,
      !    An Algorithm for Minimization Using Exact Second Derivatives,
      !    Report TP 515, Theoretical Physics Division,
      !    AERE Harwell, Oxon., England, 1973.
      !
      !    Jorge More,
      !    The Levenberg-Marquardt Algorithm, Implementation and Theory,
      !    Springer Lecture Notes in Mathematics Number 630, pages 105-116,
      !    edited by G A Watson,
      !    Springer-Verlag, Berlin and New York, 1978.
      !
      !    Jorge More and Danny Sorensen,
      !    Computing a Trust Region Step,
      !    Technical Report ANL-81-83,
      !    Argonne National Lab, 1981.
      !
      !    Richard Varga,
      !    Minimal Gerschgorin Sets,
      !    Pacific Journal of Mathematics,
      !    Volume 15, pages 719-729, 1965.
      !
      !  Parameters:
      !
      !     d (in)  = the scale vector, i.e. the diagonal of the scale
      !              matrix  d  mentioned above under purpose.
      !   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then
      !              step = 0  and  v(stppar) = 0  are returned.
      ! dihdi (in)  = lower triangle of the scaled hessian (approximation),
      !              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e.,
      !              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc.
      !    ka (i/o) = the number of hebden iterations (so far) taken to deter-
      !              mine step.  ka < 0 on input means this is the first
      !              attempt to determine step (for the present dig and dihdi)
      !              -- ka is initialized to 0 in this case.  output with
      !              ka = 0  (or v(stppar) = 0)  means  step = -(h**-1)*g.
      !     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
      !     p (in)  = number of parameters -- the hessian is a  p x p  matrix.
      !  step (i/o) = the step computed.
      !     v (i/o) contains various constants and variables described below.
      !     w (i/o) = workspace of length 4*p + 6.
      !
      !   entries in v
      !
      ! v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
      ! v(dstnrm) (output) = 2-norm of d*step.
      ! v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
      !             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
      ! v(epslon) (in)  = max. rel. error allowed for psi(step).  for the
      !             step returned, psi(step) will exceed its optimal value
      !             by less than -v(epslon)*psi(step).  suggested value = 0.1.
      ! v(gtstep) (out) = inner product between g and step.
      ! v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def.
      !             h only -- v(nreduc) is set to zero otherwise).
      ! v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
      !             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
      !             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
      ! v(phmxfc) (in)  (see v(phmnfc).)
      !             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
      ! v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
      ! v(radius) (in)  = radius of current (scaled) trust region.
      ! v(rad0)   (i/o) = value of v(radius) from previous call.
      ! v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha
      !             described below under algorithm notes.  if h + alpha*d**2
      !             (see algorithm notes) is (nearly) singular, however,
      !             then v(stppar) = -alpha.
      !
      !   usage notes
      !
      !     if it is desired to recompute step using a different value of
      !     v(radius), then this routine may be restarted by calling it
      !     with all parameters unchanged except v(radius).  (this explains
      !     why step and w are listed as i/o).  on an initial call (one with
      !     ka < 0), step and w need not be initialized and only compo-
      !     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and
      !     v(rad0) of v must be initialized.
      !
        integer ka, p
        real ( kind = 8 ) d(p), dig(p), dihdi(*), l(*), v(21), step(p), w(*)
      !     dimension dihdi(p*(p+1)/2), l(p*(p+1)/2), w(4*p+7)
      !
        logical restrt
        integer dggdmx, diag, diag0, dstsav, emax, emin, i, im1, inc, irc
        integer j, k, kalim, kamin, k1, lk0, phipin, q, q0, uk0, x
        real ( kind = 8 ) alphak, aki, akk, delta, dst, eps, gtsta, lk
        real ( kind = 8 ) oldphi
        real ( kind = 8 ) phi, phimax, phimin, psifac, rad, radsq
        real ( kind = 8 ) root, si, sk, sw, t, twopsi, t1, t2, uk, wi
        real ( kind = 8 ) big, dgxfac, epsfac, four, half, kappa, negone
        real ( kind = 8 ) one, p001, six, three, two, zero
        real ( kind = 8 ) dotprd, lsvmin, rmdcon, v2norm
        integer dgnorm, dstnrm, dst0, epslon, gtstep, stppar, nreduc
        integer phmnfc, phmxfc, preduc, radius, rad0
      
        parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19, gtstep=4 )
        parameter ( nreduc=6, phmnfc=20, phmxfc=21, preduc=7, radius=8 )
        parameter ( rad0=9, stppar=5)
      
        parameter (epsfac=50.0d+0, four=4.0d+0, half=0.5d+0 )
        parameter ( kappa=2.0d+0, negone=-1.0d+0, one=1.0d+0, p001=1.0d-3 )
        parameter ( six=6.0d+0, three=3.0d+0, two=2.0d+0, zero=0.0d+0)
      
        save dgxfac
      
        data big/0.d+0/, dgxfac/0.d+0/
      !
      !  Store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx).
      !
        dggdmx = p + 1
      !
      !  Store Gerschgorin over- and underestimates of the largest
      !  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax)
      !  and w(emin) respectively.
      !
        emax = dggdmx + 1
        emin = emax + 1
      !
      !  For use in recomputing step, the final values of lk, uk, dst,
      !  and the inverse derivative of more*s phi at 0 (for pos. def.
      !  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin)
      !  respectively.
      !
        lk0 = emin + 1
        phipin = lk0 + 1
        uk0 = phipin + 1
        dstsav = uk0 + 1
      !
      !  Store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p).
      !
        diag0 = dstsav
        diag = diag0 + 1
      !
      !  Store -d*step in w(q),...,w(q0+p).
      !
        q0 = diag0 + p
        q = q0 + 1
      !
      !  Allocate storage for scratch vector x
      !
        x = q + p
        rad = v(radius)
        radsq = rad**2
      !
      !  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of d*step.
      !
        phimax = v(phmxfc) * rad
        phimin = v(phmnfc) * rad
        psifac = two * v(epslon) / (three * (four * (v(phmnfc) + 1.0D+00 ) * &
                          (kappa + 1.0D+00 )  +  kappa  +  two) * rad**2)
      !
      !  OLDPHI is used to detect limits of numerical accuracy.  if
      !  we recompute step and it does not change, then we accept it.
      !
        oldphi = 0.0D+00
        eps = v(epslon)
        irc = 0
        restrt = .false.
        kalim = ka + 50
      !
      !  Start or restart, depending on ka
      !
        if (ka >= 0) go to 290
      !
      !  fresh start
      !
        k = 0
        uk = negone
        ka = 0
        kalim = 50
        v(dgnorm) = v2norm(p, dig)
        v(nreduc) = 0.0D+00
        v(dst0) = 0.0D+00
        kamin = 3
        if (v(dgnorm) == 0.0D+00 ) kamin = 0
      !
      !  store diag(dihdi) in w(diag0+1),...,w(diag0+p)
      !
        j = 0
        do i = 1, p
          j = j + i
          k1 = diag0 + i
          w(k1) = dihdi(j)
        end do
      !
      !  determine w(dggdmx), the largest element of dihdi
      !
        t1 = 0.0D+00
        j = p * (p + 1) / 2
        do i = 1, j
           t = abs(dihdi(i))
           if (t1 < t) t1 = t
        end do
        w(dggdmx) = t1
      !
      !  try alpha = 0
      !
       30   call lsqrt(1, p, l, dihdi, irc)
        if (irc == 0) go to 50
      !
      !  indefinite h -- underestimate smallest eigenvalue, use this
      !  estimate to initialize lower bound lk on alpha.
      !
           j = irc*(irc+1)/2
           t = l(j)
           l(j) = 1.0D+00
           w(1:irc) = 0.0D+00
           w(irc) = one
           call litvmu(irc, w, l, w)
           t1 = v2norm(irc, w)
           lk = -t / t1 / t1
           v(dst0) = -lk
           if (restrt) go to 210
           go to 70
      !
      !  positive definite h -- compute unmodified newton step.
      !
       50   lk = 0.0D+00
        t = lsvmin(p, l, w(q), w(q))
        if (t >= one) go to 60
           if (big <= 0.0D+00 ) big = rmdcon(6)
           if (v(dgnorm) >= t*t*big) go to 70
       60   call livmul(p, w(q), l, dig)
        gtsta = dotprd(p, w(q), w(q))
        v(nreduc) = half * gtsta
        call litvmu(p, w(q), l, w(q))
        dst = v2norm(p, w(q))
        v(dst0) = dst
        phi = dst - rad
        if (phi <= phimax) go to 260
        if (restrt) go to 210
      !
      !  Prepare to compute Gerschgorin estimates of largest (and
      !  smallest) eigenvalues.
      !
       70   k = 0
        do i = 1, p
           wi = 0.0D+00
           im1 = i - 1
           do j = 1, im1
             k = k + 1
             t = abs(dihdi(k))
             wi = wi + t
             w(j) = w(j) + t
           end do
           w(i) = wi
           k = k + 1
        end do
      !
      !  (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)
      !
        k = 1
        t1 = w(diag) - w(1)
      
        do i = 2, p
           j = diag0 + i
           t = w(j) - w(i)
           if ( t < t1 ) then
             t1 = t
             k = i
           end if
        end do
      
        sk = w(k)
        j = diag0 + k
        akk = w(j)
        k1 = k*(k-1)/2 + 1
        inc = 1
        t = 0.0D+00
        do i = 1, p
           if (i == k) go to 130
           aki = abs(dihdi(k1))
           si = w(i)
           j = diag0 + i
           t1 = half * (akk - w(j) + si - aki)
           t1 = t1 + sqrt(t1*t1 + sk*aki)
           if (t < t1) t = t1
           if (i < k) go to 140
       130     inc = i
       140     k1 = k1 + inc
        end do
      
        w(emin) = akk - t
        uk = v(dgnorm)/rad - w(emin)
        if (v(dgnorm) == 0.0D+00 ) uk = uk + p001 + p001*uk
        if (uk <= 0.0D+00) uk = p001
      !
      !   compute Gerschgorin overestimate of largest eigenvalue
      !
        k = 1
        t1 = w(diag) + w(1)
        if (p <= 1) go to 170
      
        do i = 2, p
           j = diag0 + i
           t = w(j) + w(i)
           if (t <= t1) go to 160
                t1 = t
                k = i
       160     continue
        end do
      
       170  sk = w(k)
        j = diag0 + k
        akk = w(j)
        k1 = k*(k-1)/2 + 1
        inc = 1
        t = 0.0D+00
      
        do i = 1, p
           if (i == k) go to 180
           aki = abs(dihdi(k1))
           si = w(i)
           j = diag0 + i
           t1 = half * (w(j) + si - aki - akk)
           t1 = t1 + sqrt(t1*t1 + sk*aki)
           if (t < t1) t = t1
           if (i < k) go to 190
       180     inc = i
       190     k1 = k1 + inc
        end do
      
        w(emax) = akk + t
        lk = max (lk, v(dgnorm)/rad - w(emax))
      !
      !  alphak = current value of alpha (see alg. notes above).  we
      !  use More's scheme for initializing it.
      !
        alphak = abs(v(stppar)) * v(rad0)/rad
      
        if (irc /= 0) go to 210
      !
      !  Compute l0 for positive definite H.
      !
        call livmul(p, w, l, w(q))
        t = v2norm(p, w)
        w(phipin) = dst / t / t
        lk = max (lk, phi*w(phipin))
      !
      !  safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)
      !
       210  ka = ka + 1
        if (-v(dst0) >= alphak .or. alphak < lk .or. alphak >= uk) then
          alphak = uk * max (p001, sqrt(lk/uk))
        end if
        if (alphak <= 0.0D+00) alphak = half * uk
        if (alphak <= 0.0D+00) alphak = uk
        k = 0
        do i = 1, p
           k = k + i
           j = diag0 + i
           dihdi(k) = w(j) + alphak
        end do
      !
      !  Try computing Cholesky decomposition
      !
        call lsqrt(1, p, l, dihdi, irc)
        if (irc == 0) go to 240
      !
      !  (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
      !  smallest eigenvalue for use in updating lk
      !
        j = (irc*(irc+1))/2
        t = l(j)
        l(j) = one
        w(1:irc) = 0.0D+00
        w(irc) = one
        call litvmu(irc, w, l, w)
        t1 = v2norm(irc, w)
        lk = alphak - t/t1/t1
        v(dst0) = -lk
        go to 210
      !
      !  Alphak makes (d**-1)*h*(d**-1) positive definite.
      !  compute q = -d*step, check for convergence.
      !
       240  call livmul(p, w(q), l, dig)
        gtsta = dotprd(p, w(q), w(q))
        call litvmu(p, w(q), l, w(q))
        dst = v2norm(p, w(q))
        phi = dst - rad
        if (phi <= phimax .and. phi >= phimin) go to 270
        if (phi == oldphi) go to 270
        oldphi = phi
        if (phi < 0.0D+00) go to 330
      !
      !  unacceptable alphak -- update lk, uk, alphak
      !
       250  if (ka >= kalim) go to 270
      !
      !  The following dmin1 is necessary because of restarts
      !
        if (phi < 0.0D+00) uk = min (uk, alphak)
      !
      !  kamin = 0 only iff the gradient vanishes
      !
        if (kamin == 0) go to 210
        call livmul(p, w, l, w(q))
        t1 = v2norm(p, w)
        alphak = alphak  +  (phi/t1) * (dst/t1) * (dst/rad)
        lk = max (lk, alphak)
        go to 210
      !
      !  Acceptable step on first try.
      !
       260  alphak = 0.0D+00
      !
      !  Successful step in general.  compute step = -(d**-1)*q
      !
       270  continue
      
        do i = 1, p
          j = q0 + i
          step(i) = -w(j)/d(i)
        end do
      
        v(gtstep) = -gtsta
        v(preduc) = half * (abs(alphak)*dst*dst + gtsta)
        go to 410
      !
      !  Restart with new radius
      !
       290  if (v(dst0) <= 0.0D+00 .or. v(dst0) - rad > phimax) go to 310
      !
      !  Prepare to return Newton step.
      !
           restrt = .true.
           ka = ka + 1
           k = 0
           do i = 1, p
             k = k + i
             j = diag0 + i
             dihdi(k) = w(j)
           end do
           uk = negone
           go to 30
      
       310  kamin = ka + 3
        if (v(dgnorm) == 0.0D+00) kamin = 0
        if (ka == 0) go to 50
      
        dst = w(dstsav)
        alphak = abs(v(stppar))
        phi = dst - rad
        t = v(dgnorm)/rad
        uk = t - w(emin)
        if (v(dgnorm) == 0.0D+00) uk = uk + p001 + p001*uk
        if (uk <= 0.0D+00) uk = p001
        if (rad > v(rad0)) go to 320
      !
      !  Smaller radius
      !
           lk = 0.0D+00
           if (alphak > 0.0D+00) lk = w(lk0)
           lk = max (lk, t - w(emax))
           if (v(dst0) > 0.0D+00) lk = max (lk, (v(dst0)-rad)*w(phipin))
           go to 250
      !
      !  Bigger radius.
      !
       320  if (alphak > 0.0D+00) uk = min (uk, w(uk0))
        lk = max (zero, -v(dst0), t - w(emax))
        if (v(dst0) > 0.0D+00) lk = max (lk, (v(dst0)-rad)*w(phipin))
        go to 250
      !
      !  Decide whether to check for special case... in practice (from
      !  the standpoint of the calling optimization code) it seems best
      !  not to check until a few iterations have failed -- hence the
      !  test on kamin below.
      !
       330  delta = alphak + min (zero, v(dst0))
        twopsi = alphak*dst*dst + gtsta
        if (ka >= kamin) go to 340
      !
      !  if the test in ref. 2 is satisfied, fall through to handle
      !  the special case (as soon as the more-sorensen test detects
      !  it).
      !
        if (delta >= psifac*twopsi) go to 370
      !
      !  check for the special case of  h + alpha*d**2  (nearly)
      !  singular.  use one step of inverse power method with start
      !  from lsvmin to obtain approximate eigenvector corresponding
      !  to smallest eigenvalue of (d**-1)*h*(d**-1).  lsvmin returns
      !  x and w with  l*w = x.
      !
       340  t = lsvmin(p, l, w(x), w)
      !
      !  Normalize w.
      !
        w(1:p) = t * w(1:p)
      !
      !  Complete current inv. power iter. -- replace w by (l**-t)*w.
      !
        call litvmu(p, w, l, w)
        t2 = one/v2norm(p, w)
      
        w(1:p) = t2 * w(1:p)
      
        t = t2 * t
      !
      !   now w is the desired approximate (unit) eigenvector and
      !   t*x = ((d**-1)*h*(d**-1) + alphak*i)*w.
      !
        sw = dotprd(p, w(q), w)
        t1 = (rad + dst) * (rad - dst)
        root = sqrt(sw*sw + t1)
        if (sw < 0.0D+00) root = -root
        si = t1 / (sw + root)
      !
      !  The actual test for the special case:
      !
        if ((t2*si)**2 <= eps*(dst**2 + alphak*radsq)) go to 380
      !
      !  Update upper bound on smallest eigenvalue (when not positive)
      !  (as recommended by more and sorensen) and continue...
      !
        if (v(dst0) <= 0.0D+00) v(dst0) = min (v(dst0), t2**2 - alphak)
        lk = max (lk, -v(dst0))
      !
      !  Check whether we can hope to detect the special case in
      !  the available arithmetic.  accept step as it is if not.
      !
      !  If not yet available, obtain machine dependent value dgxfac.
      !
       370  if (dgxfac == 0.0D+00) dgxfac = epsfac * rmdcon(3)
      
        if (delta > dgxfac*w(dggdmx)) go to 250
           go to 270
      !
      !  Special case detected... negate alphak to indicate special case
      !
       380  alphak = -alphak
        v(preduc) = half * twopsi
      !
      !  Accept current step if adding si*w would lead to a
      !  further relative reduction in psi of less than v(epslon)/3.
      !
        t1 = 0.0D+00
        t = si*(alphak*sw - half*si*(alphak + t*dotprd(p,w(x),w)))
        if (t < eps*twopsi/six) go to 390
           v(preduc) = v(preduc) + t
           dst = rad
           t1 = -si
       390  continue
      
        do i = 1, p
           j = q0 + i
           w(j) = t1*w(i) - w(j)
           step(i) = w(j) / d(i)
        end do
      
        v(gtstep) = dotprd(p, dig, w(q))
      !
      !  Save values for use in a possible restart
      !
       410  v(dstnrm) = dst
        v(stppar) = alphak
        w(lk0) = lk
        w(uk0) = uk
        v(rad0) = rad
        w(dstsav) = dst
      !
      !  Restore diagonal of DIHDI.
      !
        j = 0
        do i = 1, p
          j = j + i
          k = diag0 + i
          dihdi(j) = w(k)
        end do
      
        return
      end
      subroutine humit ( d, fx, g, h, iv, lh, liv, lv, n, v, x )
      
      !*****************************************************************************80
      !
      !! HUMIT carries out unconstrained minimization iterations for HUMSL.
      !
      !  Discussion:
      !
      !    The Hessian matrix is provided by the caller.
      !
      !    parameters iv, n, v, and x are the same as the corresponding
      !    ones to humsl (which see), except that v can be shorter (since
      !    the part of v that humsl uses for storing g and h is not needed).
      !    moreover, compared with humsl, iv(1) may have the two additional
      !    output values 1 and 2, which are explained below, as is the use
      !    of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
      !    output value from humsl, is not referenced by humit or the
      !    subroutines it calls.
      !
      !  Parameters:
      !
      ! d.... scale vector.
      ! fx... function value.
      ! g.... gradient vector.
      ! h.... lower triangle of the hessian, stored rowwise.
      ! iv... integer value array.
      ! lh... length of h = p*(p+1)/2.
      ! liv.. length of iv (at least 60).
      ! lv... length of v (at least 78 + n*(n+21)/2).
      ! n.... number of variables (components in x and g).
      ! v.... floating-point value array.
      ! x.... parameter vector.
      !
      ! iv(1) = 1 means the caller should set fx to f(x), the function value
      !             at x, and call humit again, having changed none of the
      !             other parameters.  an exception occurs if f(x) cannot be
      !             computed (e.g. if overflow would occur), which may happen
      !             because of an oversized step.  in this case the caller
      !             should set iv(toobig) = iv(2) to 1, which will cause
      !             humit to ignore fx and try a smaller step.  the para-
      !             meter nf that humsl passes to calcf (for possible use by
      !             calcgh) is a copy of iv(nfcall) = iv(6).
      ! iv(1) = 2 means the caller should set g to g(x), the gradient of f at
      !             x, and h to the lower triangle of h(x), the hessian of f
      !             at x, and call humit again, having changed none of the
      !             other parameters except perhaps the scale vector d.
      !                  the parameter nf that humsl passes to calcg is
      !             iv(nfgcal) = iv(7).  if g(x) and h(x) cannot be evaluated,
      !             then the caller may set iv(nfgcal) to 0, in which case
      !             humit will return with iv(1) = 65.
      !                  note -- humit overwrites h with the lower triangle
      !             of  diag(d)**-1 * h(x) * diag(d)**-1.
      !
        integer lh
        integer liv
        integer lv
        integer n
      
        integer iv(liv)
        real ( kind = 8 ) d(n), fx, g(n), h(lh), v(lv), x(n)
        integer dg1, i, j, k, l, lstgst, nn1o2, step1
        integer temp1, w1, x01
        real ( kind = 8 ) t
        real ( kind = 8 ) one, onep2
        logical stopx
        real ( kind = 8 ) dotprd, reldst, v2norm
        integer cnvcod, dg, dgnorm, dinit, dstnrm, dtinit, dtol
        integer dtype, d0init, f, f0, fdif, gtstep, incfac, irc, kagqt
        integer lmat, lmax0, lmaxs, mode, model, mxfcal, mxiter, nextv
        integer nfcall, nfgcal, ngcall, niter, preduc, radfac, radinc
        integer radius, rad0, reldx, restor, step, stglim, stlstg, stppar
        integer toobig, tuner4, tuner5, vneed, w, xirc, x0
      
        parameter (cnvcod=55, dg=37, dtol=59, dtype=16, irc=29, kagqt=33 )
        parameter ( lmat=42, mode=35, model=5, mxfcal=17, mxiter=18 )
        parameter ( nextv=47, nfcall=6, nfgcal=7, ngcall=30, niter=31 )
        parameter ( radinc=8, restor=9, step=40, stglim=11, stlstg=41 )
        parameter ( toobig=2, vneed=4, w=34, xirc=13, x0=43)
        parameter (dgnorm=1, dinit=38, dstnrm=2, dtinit=39, d0init=40 )
        parameter ( f=10, f0=13, fdif=11, gtstep=4, incfac=23, lmax0=35 )
        parameter ( lmaxs=36, preduc=7, radfac=16, radius=8, rad0=9 )
        parameter ( reldx=17, stppar=5, tuner4=29, tuner5=30)
      
        parameter (one=1.d+0, onep2=1.2d+0 )
      
        i = iv(1)
        if (i == 1) go to 30
        if (i == 2) go to 40
      !
      !   check validity of iv and v input values
      !
        if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
        if (iv(1) == 12 .or. iv(1) == 13) then
         iv(vneed) = iv(vneed) + n*(n+21)/2 + 7
        end if
        call parck(2, d, iv, liv, lv, n, v)
        i = iv(1) - 2
        if (i > 12) then
          return
        end if
        nn1o2 = n * (n + 1) / 2
        if (lh >= nn1o2) go to (210,210,210,210,210,210,160,120,160, &
         10,10,20),i
           iv(1) = 66
           go to 350
      !
      !   storage allocation
      !
       10   iv(dtol) = iv(lmat) + nn1o2
        iv(x0) = iv(dtol) + 2*n
        iv(step) = iv(x0) + n
        iv(stlstg) = iv(step) + n
        iv(dg) = iv(stlstg) + n
        iv(w) = iv(dg) + n
        iv(nextv) = iv(w) + 4*n + 7
        if (iv(1) /= 13) go to 20
           iv(1) = 14
           return
      !
      !   initialization
      !
       20   iv(niter) = 0
        iv(nfcall) = 1
        iv(ngcall) = 1
        iv(nfgcal) = 1
        iv(mode) = -1
        iv(model) = 1
        iv(stglim) = 1
        iv(toobig) = 0
        iv(cnvcod) = 0
        iv(radinc) = 0
        v(rad0) = 0.0D+00
        v(stppar) = 0.0D+00
        if (v(dinit) >= 0.0D+00) call vscopy(n, d, v(dinit))
        k = iv(dtol)
        if (v(dtinit) > 0.0D+00) call vscopy(n, v(k), v(dtinit))
        k = k + n
        if (v(d0init) > 0.0D+00) call vscopy(n, v(k), v(d0init))
        iv(1) = 1
        return
      
       30   v(f) = fx
        if (iv(mode) >= 0) go to 210
        iv(1) = 2
        if (iv(toobig) == 0) then
          return
        end if
           iv(1) = 63
           go to 350
      !
      !  Make sure gradient could be computed
      !
       40   if (iv(nfgcal) /= 0) go to 50
           iv(1) = 65
           go to 350
      !
      !  Update the scale vector d
      !
       50   dg1 = iv(dg)
        if (iv(dtype) <= 0) go to 70
        k = dg1
        j = 0
        do i = 1, n
           j = j + i
           v(k) = h(j)
           k = k + 1
        end do
      
        call dupdu(d, v(dg1), iv, liv, lv, n, v)
      !
      !  Compute scaled gradient and its norm
      !
       70   dg1 = iv(dg)
        k = dg1
        do i = 1, n
           v(k) = g(i) / d(i)
           k = k + 1
        end do
      
        v(dgnorm) = v2norm(n, v(dg1))
      !
      !  Compute scaled hessian
      !
        k = 1
        do i = 1, n
           t = one / d(i)
           do j = 1, i
                h(k) = t * h(k) / d(j)
                k = k + 1
           end do
        end do
      
        if (iv(cnvcod) /= 0) go to 340
        if (iv(mode) == 0) go to 300
      !
      !  Allow first step to have scaled 2-norm at most v(lmax0)
      !
        v(radius) = v(lmax0)
      
        iv(mode) = 0
      !
      !  Main loop
      !
      !  print iteration summary, check iteration limit
      !
       110  call itsum(d, g, iv, liv, lv, n, v, x)
       120  k = iv(niter)
        if (k < iv(mxiter)) go to 130
           iv(1) = 10
           go to 350
      
       130  iv(niter) = k + 1
      !
      !  initialize for start of next iteration
      !
        dg1 = iv(dg)
        x01 = iv(x0)
        v(f0) = v(f)
        iv(irc) = 4
        iv(kagqt) = -1
      !
      !  Copy x to x0
      !
        call vcopy ( n, v(x01), x )
      !
      !  Update radius
      !
        if (k == 0) go to 150
        step1 = iv(step)
        k = step1
        do i = 1, n
           v(k) = d(i) * v(k)
           k = k + 1
        end do
        v(radius) = v(radfac) * v2norm(n, v(step1))
      !
      !  Check STOPX and function evaluation limit.
      !
       150  if (.not. stopx ( ) ) go to 170
           iv(1) = 11
           go to 180
      !
      !  Come here when restarting after func. eval. limit or STOPX.
      !
       160  if (v(f) >= v(f0)) go to 170
           v(radfac) = one
           k = iv(niter)
           go to 130
      
       170  if (iv(nfcall) < iv(mxfcal)) go to 190
           iv(1) = 9
       180     if (v(f) >= v(f0)) go to 350
      !
      !  In case of STOPX or function evaluation limit with
      !  improved v(f), evaluate the gradient at x.
      !
                iv(cnvcod) = iv(1)
                go to 290
      !
      !  Compute candidate step
      !
       190  step1 = iv(step)
        dg1 = iv(dg)
        l = iv(lmat)
        w1 = iv(w)
        call gqtst(d, v(dg1), h, iv(kagqt), v(l), n, v(step1), v, v(w1))
        if (iv(irc) == 6) go to 210
      !
      !  Check whether evaluating f(x0 + step) looks worthwhile
      !
        if (v(dstnrm) <= 0.0D+00) go to 210
        if (iv(irc) /= 5) go to 200
        if (v(radfac) <= one) go to 200
        if (v(preduc) <= onep2 * v(fdif)) go to 210
      !
      !  Compute f(x0 + step)
      !
       200  x01 = iv(x0)
        step1 = iv(step)
        call vaxpy(n, x, one, v(step1), v(x01))
        iv(nfcall) = iv(nfcall) + 1
        iv(1) = 1
        iv(toobig) = 0
        return
      !
      !  Assess candidate step.
      !
       210  x01 = iv(x0)
        v(reldx) = reldst(n, d, x, v(x01))
        call assst(iv, liv, lv, v)
        step1 = iv(step)
        lstgst = iv(stlstg)
        if (iv(restor) == 1) call vcopy(n, x, v(x01))
        if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
        if (iv(restor) /= 3) go to 220
           call vcopy(n, v(step1), v(lstgst))
           call vaxpy(n, x, one, v(step1), v(x01))
           v(reldx) = reldst(n, d, x, v(x01))
      
       220  k = iv(irc)
        go to (230,260,260,260,230,240,250,250,250,250,250,250,330,300), k
      !
      !  Recompute step with new radius
      !
       230     v(radius) = v(radfac) * v(dstnrm)
           go to 150
      !
      !  Compute step of length v(lmaxs) for singular convergence test.
      !
       240  v(radius) = v(lmaxs)
        go to 190
      !
      !  Convergence or false convergence
      !
       250  iv(cnvcod) = k - 4
        if (v(f) >= v(f0)) go to 340
           if (iv(xirc) == 14) go to 340
                iv(xirc) = 14
      !
      !  Process acceptable step.
      !
       260  if (iv(irc) /= 3) go to 290
           temp1 = lstgst
      !
      !  prepare for gradient tests
      !  set  temp1 = hessian * step + g(x0)
      !  = diag(d) * (h * step + g(x0))
      !
      !  Use x0 vector as temporary.
      !
           k = x01
      
           do i = 1, n
             v(k) = d(i) * v(step1)
             k = k + 1
             step1 = step1 + 1
           end do
      
           call slvmul(n, v(temp1), h, v(x01))
      
           do i = 1, n
             v(temp1) = d(i) * v(temp1) + g(i)
             temp1 = temp1 + 1
           end do
      !
      !  Compute gradient and hessian.
      !
       290  iv(ngcall) = iv(ngcall) + 1
        iv(1) = 2
        return
      
       300  iv(1) = 2
        if (iv(irc) /= 3) go to 110
      !
      !  Set v(radfac) by gradient tests.
      !
        temp1 = iv(stlstg)
        step1 = iv(step)
      !
      !  Set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))
      !
        k = temp1
        do i = 1, n
           v(k) = ( v(k) - g(i) ) / d(i)
           k = k + 1
        end do
      !
      !  Do gradient tests,
      !
        if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) go to 320
             if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  go to 110
       320            v(radfac) = v(incfac)
                  go to 110
      !
      !  misc. details
      !
      !  bad parameters to assess
      !
       330  iv(1) = 64
        go to 350
      !
      !  Print summary of final iteration and other requested items
      !
       340  iv(1) = iv(cnvcod)
        iv(cnvcod) = 0
       350  call itsum(d, g, iv, liv, lv, n, v, x)
      
        return
      end
      subroutine humsl ( n, d, x, calcf, calcgh, iv, liv, lv, v, uiparm, &
        urparm, ufparm )
      
      !*****************************************************************************80
      !
      !! HUMSL minimizes a general unconstrained objective function.
      !
      !  Discussion:
      !
      !    The gradient and Hessian are provided by the caller.
      !
      !    this routine is like sumsl, except that the subroutine para-
      !    meter calcg of sumsl (which computes the gradient of the objec-
      !    tive function) is replaced by the subroutine parameter calcgh,
      !    which computes both the gradient and (lower triangle of the)
      !    hessian of the objective function.
      !
      !  Reference:
      !
      !    David Gay,
      !    Computing Optimal Locally Constrained Steps,
      !    SIAM Journal on Scientific and Statistical Computing,
      !    Volume 2, pages 186-197, 1981.
      !
      !  Parameters:
      !
      !    the calling sequence is...
      !             call calcgh(n, x, nf, g, h, uiparm, urparm, ufparm)
      !     parameters n, x, nf, g, uiparm, urparm, and ufparm are the same
      !     as for sumsl, while h is an array of length n*(n+1)/2 in which
      !     calcgh must store the lower triangle of the hessian at x.  start-
      !     ing at h(1), calcgh must store the hessian entries in the order
      !     (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ...
      !        the value printed (by itsum) in the column labelled stppar
      !     is the levenberg-marquardt used in computing the current step.
      !     zero means a full newton step.  if the special case described in
      !     ref. 1 is detected, then stppar is negated.  the value printed
      !     in the column labelled npreldf is zero if the current hessian
      !     is not positive definite.
      !        it sometimes proves worthwhile to let d be determined from the
      !     diagonal of the hessian matrix by setting iv(dtype) = 1 and
      !     v(dinit) = 0.  the following iv and v components are relevant...
      !
      ! iv(dtol)  iv(59) gives the starting subscript in v of the dtol
      !             array used when d is updated.  (iv(dtol) can be
      !             initialized by calling humsl with iv(1) = 13.)
      ! iv(dtype).... iv(16) tells how the scale vector d should be chosen.
      !             iv(dtype) <= 0 means that d should not be updated, and
      !             iv(dtype) >= 1 means that d should be updated as
      !             described below with v(dfac).  default = 0.
      ! v(dfac)  v(41) and the dtol and d0 arrays (see v(dtinit) and
      !             v(d0init)) are used in updating the scale vector d when
      !             iv(dtype) > 0.  (d is initialized according to
      !             v(dinit), described in sumsl.)  let
      !                  d1(i) = max(sqrt(abs(h(i,i))), v(dfac)*d(i)),
      !             where h(i,i) is the i-th diagonal element of the current
      !             hessian.  if iv(dtype) = 1, then d(i) is set to d1(i)
      !             unless d1(i) < dtol(i), in which case d(i) is set to
      !                  max(d0(i), dtol(i)).
      !             if iv(dtype) >= 2, then d is updated during the first
      !             iteration as for iv(dtype) = 1 (after any initialization
      !             due to v(dinit)) and is left unchanged thereafter.
      !             default = 0.6.
      ! v(dtinit)... v(39), if positive, is the value to which all components
      !             of the dtol array (see v(dfac)) are initialized.  if
      !             v(dtinit) = 0, then it is assumed that the caller has
      !             stored dtol in v starting at v(iv(dtol)).
      !             default = 10**-6.
      ! v(d0init)... v(40), if positive, is the value to which all components
      !             of the d0 vector (see v(dfac)) are initialized.  if
      !             v(dfac) = 0, then it is assumed that the caller has
      !             stored d0 in v starting at v(iv(dtol)+n).  default = 1.0.
      !
        integer liv
        integer lv
        integer n
      
        integer iv(liv)
        integer uiparm(*)
        real ( kind = 8 ) d(n), x(n), v(lv), urparm(*)
      !     dimension v(78 + n*(n+12)), uiparm(*), urparm(*)
        external ufparm
        integer g1, h1, iv1, lh, nf
        real ( kind = 8 ) f
        integer g, h, nextv, nfcall, nfgcal, toobig, vneed
      
        parameter (nextv=47, nfcall=6, nfgcal=7, g=28, h=56, toobig=2, &
      vneed=4)
      
        lh = n * (n + 1) / 2
      
        if ( iv(1) == 0 ) then
          call deflt ( 2, iv, liv, lv, v )
        end if
      
        if (iv(1) == 12 .or. iv(1) == 13) then
          iv(vneed) = iv(vneed) + n*(n+3)/2
        end if
      
        iv1 = iv(1)
        if (iv1 == 14) go to 10
        if (iv1 > 2 .and. iv1 < 12) go to 10
        g1 = 1
        h1 = 1
        if (iv1 == 12) iv(1) = 13
        go to 20
      
       10   g1 = iv(g)
        h1 = iv(h)
      
       20   call humit(d, f, v(g1), v(h1), iv, lh, liv, lv, n, v, x)
      
        if (iv(1) - 2) 30, 40, 50
      
       30   nf = iv(nfcall)
        call calcf(n, x, nf, f, uiparm, urparm, ufparm)
        if (nf <= 0) iv(toobig) = 1
        go to 20
      
       40   call calcgh(n, x, iv(nfgcal), v(g1), v(h1), uiparm, &
       urparm, ufparm)
        go to 20
      
       50   if (iv(1) /= 14) then
              return
            end if
      !
      !  storage allocation
      !
        iv(g) = iv(nextv)
        iv(h) = iv(g) + n
        iv(nextv) = iv(h) + n*(n+1)/2
      
        if ( iv1 /= 13 ) then
          go to 10
        end if
      
        return
      end
      subroutine itsum ( d, g, iv, liv, lv, p, v, x )
  
      !use m_optimization
      !
      !*****************************************************************************80
      !
      !! ITSUM prints an iteration summary.
      !
        integer liv
        integer lv
        integer p
      
        real ( kind = 8 ) d(p)
        real ( kind = 8 ) g(p)
        integer iv(liv)
        real ( kind = 8 ) v(lv)
        real ( kind = 8 ) x(p)
        integer alg, i, iv1, m, nf, ng, ol, pu
        character*4 model1(6), model2(6)
        real ( kind = 8 ) nreldf, oldf, preldf, reldf
        integer algsav, dstnrm, f, fdif, f0, needhd, nfcall, nfcov, ngcov
        integer ngcall, niter, nreduc, outlev, preduc, prntit, prunit
        integer reldx, solprt, statpr, stppar, sused, x0prt
      
        parameter (algsav=51, needhd=36, nfcall=6, nfcov=52, ngcall=30 )
        parameter ( ngcov=53, niter=31, outlev=19, prntit=39, prunit=21 )
        parameter ( solprt=22, statpr=23, sused=64, x0prt=24)
        parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6, preduc=7, &
      reldx=17 )
        parameter ( stppar=5)
      
        data model1/'    ','    ','    ','    ','  g ','  s '/
        data model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/
      
        pu = iv(prunit)
      
        if ( pu == 0 ) then
          return
        end if
      
        iv1 = iv(1)
        if (iv1 > 62) iv1 = iv1 - 51
        ol = iv(outlev)
        alg = iv(algsav)
        if (iv1 < 2 .or. iv1 > 15) go to 370
        if (iv1 >= 12) go to 120
        if (iv1 == 2 .and. iv(niter) == 0) go to 390
        if (ol == 0) go to 120
        if (iv1 >= 10 .and. iv(prntit) == 0) go to 120
        if (iv1 > 2) go to 10
           iv(prntit) = iv(prntit) + 1
           if (iv(prntit) < iabs(ol)) then
             return
           end if
   10   nf = iv(nfcall) - iabs(iv(nfcov))
        iv(prntit) = 0
        reldf = 0.0D+00
        preldf = 0.0D+00
        oldf = max (abs(v(f0)), abs(v(f)))
        if (oldf <= 0.0D+00) go to 20
           reldf = v(fdif) / oldf
           !Andy:
       !   print *,'v(fdif),oldf='
       !   write(*,'(2F25.20)') v(fdif),oldf,v(f),v(f0)
           preldf = v(preduc) / oldf
   20   if (ol > 0) go to 60
      !
      !  print short summary line
      !
           if (iv(needhd) == 1 .and. alg == 1) write(pu,30)
   30   format(/'   it   nf',6x,'f',7x,'reldf',3x,'preldf', & 
           3x,'reldx',2x,'model  stppar')
           if (iv(needhd) == 1 .and. alg == 2) write(pu,40)
   40  format(/'    it   nf',7x,'f',8x,'reldf',4x,'preldf',4x, &
           'reldx  stppar')
           iv(needhd) = 0
           if (alg == 2) go to 50
           m = iv(sused)
           write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
             model1(m), model2(m), v(stppar)
       !   !Andy: Check if optfunc is less than ... after ... fn evals.
       !   if ((nf.gt.n_optdiffint).and.(v(f).gt.optdiffint)) then
             ! write(*,'(A9,F15.10,A14,F15.10,A7,I5,A26)') &
             !       ' Optfunc=',v(f),' greater than ',optdiffint, &
             !       ' after ',nf,' function evals: Stopping.'
       !      opt_failed=.true.
       !      return
       !   endif
           go to 120
      
   50  write(pu,110) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
      v(stppar)
           go to 120
      !
      !  print long summary line
      !
   60   if (iv(needhd) == 1 .and. alg == 1) write(pu,70)
   70   format(/11h    it   nf,6x,1hf,7x,5hreldf,3x,6hpreldf,3x,5hreldx, &
              2x,13hmodel  stppar,2x,6hd*step,2x,7hnpreldf)
        if (iv(needhd) == 1 .and. alg == 2) write(pu,80)
   80   format(/11h    it   nf,7x,1hf,8x,5hreldf,4x,6hpreldf,4x,5hreldx, &
              3x,6hstppar,3x,6hd*step,3x,7hnpreldf)
        iv(needhd) = 0
        nreldf = 0.0D+00
        if (oldf > 0.0D+00) nreldf = v(nreduc) / oldf
        if (alg == 2) go to 90
        m = iv(sused)
        write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
                    model1(m), model2(m), v(stppar), v(dstnrm), nreldf
        go to 120
      
   90   write(pu,110) iv(niter), nf, v(f), reldf, preldf, &
               v(reldx), v(stppar), v(dstnrm), nreldf
       !Andy: Check if optfunc is less than ... after ... fn evals.
      !if ((nf.gt.n_optdiffint).and.(v(f).gt.optdiffint)) then
      !   write(*,'(A9,F15.10,A14,F15.10,A7,I5,A26)') &
      !         ' Optfunc=',v(f),' greater than ',optdiffint, &
      !         ' after ',nf,' function evals: Stopping.'
      !   opt_failed=.true.
      !   return
      !endif
  100  format(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)
!  110  format(i6,i5,d25.15,2d25.15,3d25.15,d25.15)
  110  format(i6,i5,d11.3,2d10.2,3d9.1,d10.2)
     
  120  if (iv(statpr) < 0) go to 430
        go to (999, 999, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, &
          330, 350, 520), iv1
      
  130  write(pu,140)
  140  format(/' x-convergence')
        go to 430
      
  150  write ( pu, 160 )
  160  format(/' relative function convergence')
        go to 430
      
  170  write(pu,180)
  180  format(/' x- and relative function convergence')
        go to 430
      
  190  write(pu,200)
  200  format(/' Absolute function convergence.')
        go to 430
      
  210  write(pu,220)
  220  format(/' Singular convergence.')
        go to 430
      
  230  write(pu,240)
  240  format(/' Convergence within error of optimization function.')
        go to 430
      
  250  write(pu,260)
  260  format(/' Function evaluation limit.')
        go to 430
      
  270  write(pu,280)
  280  format(/' Iteration limit.')
        go to 430
      
  290  write(pu,300)
  300  format(/' STOPX')
        go to 430
      
  310  write(pu,320)
  320  format(/' Initial f(x) cannot be computed.')
      
        go to 390
      
  330  write(pu,340)
  340  format(/' Bad parameters to assess.')
        go to 999
      
  350  write(pu,360)
  360  format(/' Gradient could not be computed.')
        if (iv(niter) > 0) go to 480
        go to 390
      
  370  write(pu,380) iv(1)
  380  format(/' iv(1) =',i5)
        go to 999
      !
      !   initial call on itsum
      !
  390  if (iv(x0prt) /= 0) write(pu,400) (i, x(i), d(i), i = 1, p)
  400  format(/23h     i     initial x(i),8x,4hd(i)//(1x,i5,d17.6,d14.3))
      !     the following are to avoid undefined variables when the
      !     function evaluation limit is 1...
      !
        v(dstnrm) = 0.0D+00
        v(fdif) = 0.0D+00
        v(nreduc) = 0.0D+00
        v(preduc) = 0.0D+00
        v(reldx)  = 0.0D+00
        if (iv1 >= 12) go to 999
        iv(needhd) = 0
        iv(prntit) = 0
        if (ol == 0) go to 999
        if (ol < 0 .and. alg == 1) write(pu,30)
        if (ol < 0 .and. alg == 2) write(pu,40)
        if (ol > 0 .and. alg == 1) write(pu,70)
        if (ol > 0 .and. alg == 2) write(pu,80)
        if (alg == 1) write(pu,410) v(f)
        if (alg == 2) write(pu,420) v(f)
  410  format(/11h     0    1,d10.3)
      !365  format(/11h     0    1,e11.3)
  420  format(/11h     0    1,d11.3)
        go to 999
      !
      !  Print various information requested on solution
      !
  430  iv(needhd) = 1
        if (iv(statpr) == 0) go to 480
           oldf = max (abs(v(f0)), abs(v(f)))
           preldf = 0.0D+00
           nreldf = 0.0D+00
           if (oldf <= 0.0D+00) go to 440
                preldf = v(preduc) / oldf
                nreldf = v(nreduc) / oldf
  440     nf = iv(nfcall) - iv(nfcov)
           ng = iv(ngcall) - iv(ngcov)
           write(pu,450) v(f), v(reldx), nf, ng, preldf, nreldf
  450  format(/9h function,d17.6,8h   reldx,d17.3/12h func. evals, &
          i8,9x,11hgrad. evals,i8/7h preldf,d16.3,6x,7hnpreldf,d15.3)
      
           if (iv(nfcov) > 0) write(pu,460) iv(nfcov)
  460     format(/1x,i4,50h extra func. evals for covariance and &
      diagnostics.)
           if (iv(ngcov) > 0) write(pu,470) iv(ngcov)
  470     format(1x,i4,50h extra grad. evals for covariance and &
      diagnostics.)
      
  480  if (iv(solprt) == 0) go to 999
           iv(needhd) = 1
           write(pu,490)
  490  format(/22h     i      final x(i),8x,4hd(i),10x,4hg(i)/)
           do i = 1, p
                write(pu,510) i, x(i), d(i), g(i)
           end do
  510     format(1x,i5,d16.6,2d14.3)
        go to 999
      
  520  write(pu,530)
  530  format(/'Inconsistent dimensions.')
  999  continue
      
        return
      end
      subroutine litvmu ( n, x, l, y )
      
      !*****************************************************************************80
      !
      !! LITVMU solves L' * x = y.
      !
      !  Discussion:
      !
      !    L is an  n x n  lower triangular
      !    matrix stored compactly by rows.  x and y may occupy the same
      !    storage.
      !
        integer n
      
        real ( kind = 8 ) l(*)
        real ( kind = 8 ) x(n)
        real ( kind = 8 ) y(n)
        integer i, ii, ij, i0, j
        real ( kind = 8 ) xi
      
        x(1:n) = y(1:n)
      
        i0 = n*(n+1)/2
      
        do ii = 1, n
          i = n+1 - ii
          xi = x(i)/l(i0)
          x(i) = xi
          if ( i <= 1 ) then
            exit
          end if
          i0 = i0 - i
          if ( xi /= 0.0D+00 ) then
            do j = 1, i-1
              ij = i0 + j
              x(j) = x(j) - xi*l(ij)
            end do
          end if
        end do
      
        return
      end
      subroutine livmul ( n, x, l, y )
      
      !*****************************************************************************80
      !
      !! LIVMUL solves L * x = y.
      !
      !  Discussion:
      !
      !    L is an  n x n  lower triangular
      !    matrix stored compactly by rows.  x and y may occupy the same
      !    storage.
      !
        integer n
      
        real ( kind = 8 ) x(n), l(*), y(n)
        external dotprd
        real ( kind = 8 ) dotprd
        integer i, j, k
        real ( kind = 8 ) t
      
        do k = 1, n
          if (y(k) /= 0.0D+00 ) go to 20
          x(k) = 0.0D+00
        end do
      
        return
      
  20 continue
      
        j = k*(k+1)/2
        x(k) = y(k) / l(j)
      
        if (k >= n) then
          return
        end if
      
        k = k + 1
      
        do i = k, n
           t = dotprd(i-1, l(j+1), x)
           j = j + i
           x(i) = (y(i) - t)/l(j)
        end do
      
        return
      end
      subroutine lsqrt ( n1, n, l, a, irc )
      
      !*****************************************************************************80
      !
      !! LSQRT computes rows N1 through N of the Cholesky factor L.
      !
      !  Discussion:
      !
      !   The Cholesky factor L satisfies a = l*(l**t),  where  l  and the
      !   lower triangle of  a  are both stored compactly by rows (and may
      !   occupy
      !   the same storage).
      !
      !   irc = 0 means all went well.  irc = j means the leading
      !   principal  j x j  submatrix of  a  is not positive definite --
      !   and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal.
      !
      !  Parameters:
      !
        integer n1
        integer n, irc
        real ( kind = 8 ) l(*), a(*)
      !     dimension l(n*(n+1)/2), a(n*(n+1)/2)
      !
        integer i, ij, ik, im1, i0, j, jk, j0, k
        real ( kind = 8 ) t, td
      
        i0 = n1 * (n1 - 1) / 2
      
        do i = n1, n
           td = 0.0D+00
           if (i == 1) go to 40
           j0 = 0
           im1 = i - 1
           do j = 1, im1
                t = 0.0D+00
                do k = 1, j-1
                  ik = i0 + k
                  jk = j0 + k
                  t = t + l(ik)*l(jk)
                end do
                ij = i0 + j
                j0 = j0 + j
                t = (a(ij) - t) / l(j0)
                l(ij) = t
                td = td + t*t
             end do
   40    i0 = i0 + i
           t = a(i0) - td
           if (t <= 0.0D+00) go to 60
           l(i0) = sqrt(t)
        end do
      
        irc = 0
        return
      
   60   l(i0) = t
        irc = i
      
        return
      end
      function lsvmin ( p, l, x, y )
      
      !*****************************************************************************80
      !
      !! LSVMIN estimates the smallest singular value of matrix L.
      !
      !  Discussion:
      !
      !    L is a packed lower triangular matrix.
      !
      !    this function returns a good overestimate of the smallest
      !    singular value of the packed lower triangular matrix l.
      !
      !  Reference:
      !
      !    Alan Cline, Cleve Moler, G Stewart, and James Wilkinson,
      !    An estimate for the Condition Number of a Matrix,
      !    Report TM-310, 1977,
      !    Applied Mathematics Division,
      !    Argonne National Laboratory.
      !
      !    D C Hoaglin,
      !    Theoretical properties of congruential random-number generators,
      !    an empirical view,
      !    memorandum ns-340, dept. of statistics, harvard univ., 1976.
      !
      !    D E Knuth,
      !    The Art of Computer Programming,
      !    Volume 2: Seminumerical Algorithms,
      !    Addison-wesley, reading, mass., 1969.
      !
      !    C S Smith,
      !    Multiplicative pseudo-random number generators with prime modulus,
      !    Journal of the Association for Computing Machinery,
      !    Volume 18, pages 586-593, 1971.
      !
      !  Parameters:
      !
      !  p (in)  = the order of l.  l is a  p x p  lower triangular matrix.
      !  l (in)  = array holding the elements of  l  in row order, i.e.
      !             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc.
      !  x (out) if lsvmin returns a positive value, then x is a normalized
      !             approximate left singular vector corresponding to the
      !             smallest singular value.  this approximation may be very
      !             crude.  if lsvmin returns zero, then some components of x
      !             are zero and the rest retain their input values.
      !  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an
      !             unnormalized approximate right singular vector correspond-
      !             ing to the smallest singular value.  this approximation
      !             may be crude.  if lsvmin returns zero, then y retains its
      !             input value.  the caller may pass the same vector for x
      !             and y (nonstandard fortran usage), in which case y over-
      !             writes x (for nonzero lsvmin returns).
      !
      !   algorithm notes
      !
      !     the algorithm is based on (1), with the additional provision that
      !     lsvmin = 0 is returned if the smallest diagonal element of l
      !     (in magnitude) is not more than the unit roundoff times the
      !     largest.  the algorithm uses a random number generator proposed
      !     in (4), which passes the spectral test with flying colors -- see
      !     (2) and (3).
      !
        integer p
        real ( kind = 8 ) lsvmin
        real ( kind = 8 ) l(*), x(p), y(p)
      !     dimension l(p*(p+1)/2)
      !
        integer i, ii, ix, j, ji, jj, jjj, jm1, j0, pm1
        real ( kind = 8 ) b, sminus, splus, t, xminus, xplus
        real ( kind = 8 ) half, one, r9973
        real ( kind = 8 ) dotprd, v2norm
      
        parameter (half=0.5d+0, one=1.d+0, r9973=9973.d+0 )
      
        ix = 2
        pm1 = p - 1
      !
      !  First check whether to return lsvmin = 0 and initialize x
      !
        ii = 0
        j0 = p*pm1/2
        jj = j0 + p
        if (l(jj) == 0.0D+00) go to 110
        ix = mod(3432*ix, 9973)
        b = half*(one + float(ix)/r9973)
        xplus = b / l(jj)
        x(p) = xplus
        if (p <= 1) go to 60
        do i = 1, pm1
           ii = ii + i
           if (l(ii) == 0.0D+00) go to 110
           ji = j0 + i
           x(i) = xplus * l(ji)
        end do
      !
      !  Solve (l**t)*x = b, where the components of b have randomly
      !  chosen magnitudes in (.5,1) with signs chosen to make x large.
      !
      !     do j = p-1 to 1 by -1...
        do 50 jjj = 1, pm1
           j = p - jjj
      !
      !  determine x(j) in this iteration. note for i = 1,2,...,j
      !  that x(i) holds the current partial sum for row i.
      !
           ix = mod(3432*ix, 9973)
           b = half*(one + float(ix)/r9973)
           xplus = (b - x(j))
           xminus = (-b - x(j))
           splus = abs(xplus)
           sminus = abs(xminus)
           jm1 = j - 1
           j0 = j*jm1/2
           jj = j0 + j
           xplus = xplus/l(jj)
           xminus = xminus/l(jj)
      
           do i = 1, jm1
                ji = j0 + i
                splus = splus + abs(x(i) + l(ji)*xplus)
                sminus = sminus + abs(x(i) + l(ji)*xminus)
           end do
      
   30      if (sminus > splus) xplus = xminus
           x(j) = xplus
      !
      !  update partial sums
      !
           if (jm1 > 0) call vaxpy(jm1, x, xplus, l(j0+1), x)
   50      continue
      !
      !  normalize x
      !
   60   t = one/v2norm(p, x)
      
        x(1:p) = t * x(1:p)
      !
      !  solve l*y = x and return lsvmin = 1/twonorm(y)
      !
        do j = 1, p
           jm1 = j - 1
           j0 = j*jm1/2
           jj = j0 + j
           t = 0.0D+00
           if (jm1 > 0) t = dotprd(jm1, l(j0+1), y)
           y(j) = (x(j) - t) / l(jj)
        end do
      
        lsvmin = one/v2norm(p, y)
        return
      
  110  lsvmin = 0.0D+00
        return
      end
      subroutine ltvmul ( n, x, l, y )
      
      !*****************************************************************************80
      !
      !! LTVMUL computes  x = (l**t)*y.
      !
      !  Discussion:
      !
      !    L is an  n x n  lower triangular matrix stored compactly by rows.
      !    x and y may occupy the same storage.
      !
        integer n
        real ( kind = 8 ) x(n), l(*), y(n)
      !     dimension l(n*(n+1)/2)
        integer i, ij, i0, j
        real ( kind = 8 ) yi
      
        i0 = 0
        do i = 1, n
          yi = y(i)
          x(i) = 0.0D+00
          do j = 1, i
            ij = i0 + j
            x(j) = x(j) + yi * l(ij)
          end do
          i0 = i0 + i
        end do
      
        return
      end
      subroutine lupdat ( beta, gamma, l, lambda, lplus, n, w, z )
      
      !*****************************************************************************80
      !
      !! LUPDAT computes lplus = secant update of L.
      !
      !  Discussion:
      !
      !    this routine updates the cholesky factor  l  of a symmetric
      !    positive definite matrix to which a secant update is being
      !    applied -- it computes a cholesky factor  lplus  of
      !    l * (i + z*w**t) * (i + w*z**t) * l**t.  it is assumed that  w
      !    and  z  have been chosen so that the updated matrix is strictly
      !    positive definite.
      !
      !    this code uses recurrence 3 of ref. 1 (with d(j) = 1 for all j)
      !    to compute  lplus  of the form  l * (i + z*w**t) * q,  where  q
      !    is an orthogonal matrix that makes the result lower triangular.
      !    lplus may have some negative diagonal elements.
      !
      !  Reference:
      !
      !    D Goldfarb,
      !    Factorized Variable Metric Methods for Unconstrained Optimization,
      !    Mathematics of Computation,
      !    Volume 30, pages 796-811, 1976.
      !
      !  Parameters:
      !
      !   beta = scratch vector.
      !  gamma = scratch vector.
      !      l (input) lower triangular matrix, stored rowwise.
      ! lambda = scratch vector.
      !  lplus (output) lower triangular matrix, stored rowwise, which may
      !             occupy the same storage as  l.
      !      n (input) length of vector parameters and order of matrices.
      !      w (input, destroyed on output) right singular vector of rank 1
      !             correction to  l.
      !      z (input, destroyed on output) left singular vector of rank 1
      !             correction to  l.
      !
        integer n
        real ( kind = 8 ) beta(n), gamma(n), l(*), lambda(n), lplus(*), &
         w(n),z(n)
      !     dimension l(n*(n+1)/2), lplus(n*(n+1)/2)
      !
        integer i, ij, j, jj, jp1, k, nm1
        integer np1
        real ( kind = 8 ) a, b, bj, eta, gj, lj, lij, ljj, nu, s, &
         theta, wj,zj
        real ( kind = 8 ) one
      
        parameter (one=1.d+0 )
      
        nu = one
        eta = 0.0D+00
        if (n <= 1) go to 30
        nm1 = n - 1
      !
      !  temporarily store s(j) = sum over k = j+1 to n of w(k)**2 in
      !  lambda(j).
      !
        s = 0.0D+00
        do i = 1, nm1
           j = n - i
           s = s + w(j+1)**2
           lambda(j) = s
        end do
      !
      !  compute lambda, gamma, and beta by goldfarb*s recurrence 3.
      !
        do 20 j = 1, nm1
           wj = w(j)
           a = nu*z(j) - eta*wj
           theta = one + a*wj
           s = a*lambda(j)
           lj = sqrt(theta**2 + a*s)
           if (theta > 0.0D+00) lj = -lj
           lambda(j) = lj
           b = theta*wj + s
           gamma(j) = b * nu / lj
           beta(j) = (a - b*eta) / lj
           nu = -nu / lj
           eta = -(eta + (a**2)/(theta - lj)) / lj
  20      continue
  30   lambda(n) = one + (nu*z(n) - eta*w(n))*w(n)
      !
      !  update l, gradually overwriting  w  and  z  with  l*w  and  l*z.
      !
        np1 = n + 1
        jj = n * (n + 1) / 2
      
        do k = 1, n
      
           j = np1 - k
           lj = lambda(j)
           ljj = l(jj)
           lplus(jj) = lj * ljj
           wj = w(j)
           w(j) = ljj * wj
           zj = z(j)
           z(j) = ljj * zj
           if (k == 1) go to 50
           bj = beta(j)
           gj = gamma(j)
           ij = jj + j
           jp1 = j + 1
      
           do i = jp1, n
                lij = l(ij)
                lplus(ij) = lj*lij + bj*w(i) + gj*z(i)
                w(i) = w(i) + lij*wj
                z(i) = z(i) + lij*zj
                ij = ij + i
           end do
      
   50      jj = jj - j
      
        end do
      
        return
      end
      subroutine lvmul ( n, x, l, y )
      
      !*****************************************************************************80
      !
      !! LVMUL computes x = L * y.
      !
      !  Discussion:
      !
      !    L  is an  n x n  lower triangular matrix stored compactly by rows.
      !    x and y may occupy the same storage.
      !
        integer n
      
        real ( kind = 8 ) x(n), l(*), y(n)
      !     dimension l(n*(n+1)/2)
        integer i, ii, ij, i0, j, np1
        real ( kind = 8 ) t
      
        np1 = n + 1
        i0 = n*(n+1)/2
      
        do ii = 1, n
          i = np1 - ii
          i0 = i0 - i
          t = 0.0D+00
          do j = 1, i
            ij = i0 + j
            t = t + l(ij)*y(j)
          end do
          x(i) = t
        end do
      
        return
      end
      subroutine parck ( alg, d, iv, liv, lv, n, v )
      
      !*****************************************************************************80
      !
      !! PARCK checks parameters, prints changed values.
      !
      !  Discussion:
      !
      !    alg = 1 for regression, alg = 2 for general unconstrained opt.
      !
        integer alg, liv, lv, n
        integer iv(liv)
        real ( kind = 8 ) d(n), v(lv)
        real ( kind = 8 ) rmdcon
        integer max0
        integer i, ii, iv1, j, k, l, m, miv1, miv2, ndfalt, parsv1, pu
        integer ijmp, jlim(2), miniv(2), ndflt(2)
        character*1 varnm(2), sh(2)
        character*4 cngd(3), dflt(3), vn(2,34), which(3)
        real ( kind = 8 ) big, machep, tiny, vk, vm(34), vx(34)
        integer algsav, dinit, dtype, dtype0, epslon, inits, ivneed
        integer lastiv, lastv, lmat, nextiv, nextv, nvdflt, oldn
        integer parprt, parsav, perm, prunit, vneed
      
        parameter (algsav=51, dinit=38, dtype=16, dtype0=54, epslon=19 )
        parameter ( inits=25, ivneed=3, lastiv=44, lastv=45, lmat=42 )
        parameter ( nextiv=46, nextv=47, nvdflt=50, oldn=38, parprt=20 )
        parameter ( parsav=49, perm=58, prunit=21, vneed=4)
        save big, machep, tiny
      
        data big/0.d+0/, machep/-1.d+0/, tiny/1.d+0/
      
           data vn(1,1),vn(2,1)/'epsl','on..'/
           data vn(1,2),vn(2,2)/'phmn','fc..'/
           data vn(1,3),vn(2,3)/'phmx','fc..'/
           data vn(1,4),vn(2,4)/'decf','ac..'/
           data vn(1,5),vn(2,5)/'incf','ac..'/
           data vn(1,6),vn(2,6)/'rdfc','mn..'/
           data vn(1,7),vn(2,7)/'rdfc','mx..'/
           data vn(1,8),vn(2,8)/'tune','r1..'/
           data vn(1,9),vn(2,9)/'tune','r2..'/
           data vn(1,10),vn(2,10)/'tune','r3..'/
           data vn(1,11),vn(2,11)/'tune','r4..'/
           data vn(1,12),vn(2,12)/'tune','r5..'/
           data vn(1,13),vn(2,13)/'afct','ol..'/
           data vn(1,14),vn(2,14)/'rfct','ol..'/
           data vn(1,15),vn(2,15)/'xcto','l...'/
           data vn(1,16),vn(2,16)/'xfto','l...'/
           data vn(1,17),vn(2,17)/'lmax','0...'/
           data vn(1,18),vn(2,18)/'lmax','s...'/
           data vn(1,19),vn(2,19)/'scto','l...'/
           data vn(1,20),vn(2,20)/'dini','t...'/
           data vn(1,21),vn(2,21)/'dtin','it..'/
           data vn(1,22),vn(2,22)/'d0in','it..'/
           data vn(1,23),vn(2,23)/'dfac','....'/
           data vn(1,24),vn(2,24)/'dltf','dc..'/
           data vn(1,25),vn(2,25)/'dltf','dj..'/
           data vn(1,26),vn(2,26)/'delt','a0..'/
           data vn(1,27),vn(2,27)/'fuzz','....'/
           data vn(1,28),vn(2,28)/'rlim','it..'/
           data vn(1,29),vn(2,29)/'cosm','in..'/
           data vn(1,30),vn(2,30)/'hube','rc..'/
           data vn(1,31),vn(2,31)/'rspt','ol..'/
           data vn(1,32),vn(2,32)/'sigm','in..'/
           data vn(1,33),vn(2,33)/'eta0','....'/
           data vn(1,34),vn(2,34)/'bias','....'/
      
        data vm(1)/1.0d-3/, vm(2)/-0.99d+0/, vm(3)/1.0d-3/, &
             vm(4)/1.0d-2/
        data vm(5)/1.2d+0/, vm(6)/1.d-2/, vm(7)/1.2d+0/, vm(8)/0.d+0/
        data vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(13)/0.d+0/
        data vm(15)/0.d+0/, vm(16)/0.d+0/, vm(19)/0.d+0/, &
             vm(20)/-10.d+0/
        data vm(21)/0.d+0/, vm(22)/0.d+0/, vm(23)/0.d+0/, &
             vm(27)/1.01d+0/
        data vm(28)/1.d+10/, vm(30)/0.d+0/, vm(31)/0.d+0/, vm(32)/0.d+0/
        data vm(34)/0.d+0/
      
        data vx(1)/0.9d+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8d+0/
        data vx(5)/1.d+2/, vx(6)/0.8d+0/, vx(7)/1.d+2/, vx(8)/0.5d+0/
        data vx(9)/0.5d+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1d+0/
        data vx(15)/1.d+0/, vx(16)/1.d+0/, vx(19)/1.d+0/, vx(23)/1.d+0/
        data vx(24)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+10/
        data vx(29)/1.d+0/, vx(31)/1.d+0/, vx(32)/1.d+0/, vx(33)/1.d+0/
        data vx(34)/1.d+0/
      
        data varnm(1)/'p'/, varnm(2)/'n'/, sh(1)/'s'/, sh(2)/'h'/
        data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/
        data dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
        data ijmp/33/, jlim(1)/0/, jlim(2)/24/, ndflt(1)/32/, ndflt(2)/25/
        data miniv(1)/80/, miniv(2)/59/
      
        pu = 0
        if (prunit <= liv) pu = iv(prunit)
        if (alg < 1 .or. alg > 2) go to 340
        if (iv(1) == 0) call deflt(alg, iv, liv, lv, v)
        iv1 = iv(1)
        if (iv1 /= 13 .and. iv1 /= 12) go to 10
        miv1 = miniv(alg)
        if (perm <= liv) miv1 = max0(miv1, iv(perm) - 1)
        if (ivneed <= liv) miv2 = miv1 + max0(iv(ivneed), 0)
        if (lastiv <= liv) iv(lastiv) = miv2
        if (liv < miv1) go to 300
        iv(ivneed) = 0
        iv(lastv) = max0(iv(vneed), 0) + iv(lmat) - 1
        iv(vneed) = 0
        if (liv < miv2) go to 300
        if (lv < iv(lastv)) go to 320
   10   if (alg == iv(algsav)) go to 30
        if (pu /= 0) write(pu,20) alg, iv(algsav)
   20 format(/39h the first parameter to deflt should be,i3,12h &
            rather than,i3)
           iv(1) = 82
           return
   30   if (iv1 < 12 .or. iv1 > 14) go to 60
           if (n >= 1) go to 50
                iv(1) = 81
                if (pu == 0) then
                  return
                end if
                write(pu,40) varnm(alg), n
   40           format(/8h /// bad,a1,2h =,i5)
                return
   50      if (iv1 /= 14) iv(nextiv) = iv(perm)
           if (iv1 /= 14) iv(nextv) = iv(lmat)
           if (iv1 == 13) then
             return
           end if
           k = iv(parsav) - epslon
           call vdflt(alg, lv-k, v(k+1))
           iv(dtype0) = 2 - alg
           iv(oldn) = n
           which(1) = dflt(1)
           which(2) = dflt(2)
           which(3) = dflt(3)
           go to 110
   60   if (n == iv(oldn)) go to 80
           iv(1) = 17
           if (pu == 0) then
             return
           end if
           write(pu,70) varnm(alg), iv(oldn), n
   70      format(/5h /// ,1a1,14h changed from ,i5,4h to ,i5)
           return
      
   80   if (iv1 <= 11 .and. iv1 >= 1) go to 100
           iv(1) = 80
           if (pu /= 0) write(pu,90) iv1
   90      format(/13h ///  iv(1) =,i5,28h should be between 0 and 14.)
           return
      
  100  which(1) = cngd(1)
        which(2) = cngd(2)
        which(3) = cngd(3)
      
  110  if (iv1 == 14) iv1 = 12
        if (big > tiny) go to 120
           tiny = rmdcon(1)
           machep = rmdcon(3)
           big = rmdcon(6)
           vm(12) = machep
           vx(12) = big
           vx(13) = big
           vm(14) = machep
           vm(17) = tiny
           vx(17) = big
           vm(18) = tiny
           vx(18) = big
           vx(20) = big
           vx(21) = big
           vx(22) = big
           vm(24) = machep
           vm(25) = machep
           vm(26) = machep
           vx(28) = rmdcon(5)
           vm(29) = machep
           vx(30) = big
           vm(33) = machep
  120  m = 0
        i = 1
        j = jlim(alg)
        k = epslon
        ndfalt = ndflt(alg)
      
        do l = 1, ndfalt
          vk = v(k)
          if (vk >= vm(i) .and. vk <= vx(i)) go to 140
            m = k
            if (pu /= 0) write(pu,130) vn(1,i), vn(2,i), k, vk,vm(i), &
         vx(i)
  130  format(/6h ///  ,2a4,5h.. v(,i2,3h) =,d11.3,7h should, &
            11h be between,d11.3,4h and,d11.3)
  140  k = k + 1
           i = i + 1
           if (i == j) i = ijmp
        end do
      
        if (iv(nvdflt) == ndfalt) go to 170
           iv(1) = 51
           if (pu == 0) then
             return
           end if
           write(pu,160) iv(nvdflt), ndfalt
  160     format(/13h iv(nvdflt) =,i5,13h rather than ,i5)
           return
  170  if ((iv(dtype) > 0 .or. v(dinit) > 0.0D+00) .and. iv1 == 12) then
                   go to 200
        end if
      
        do i = 1, n
           if (d(i) > 0.0D+00) go to 190
                m = 18
                if (pu /= 0) write(pu,180) i, d(i)
  180     format(/8h ///  d(,i3,3h) =,d11.3,19h should be positive)
  190     continue
        end do
      
  200  if (m == 0) go to 210
           iv(1) = m
           return
      
  210  if (pu == 0 .or. iv(parprt) == 0) then
              return
            end if
        if (iv1 /= 12 .or. iv(inits) == alg-1) go to 230
           m = 1
           write(pu,220) sh(alg), iv(inits)
  220 format(/22h nondefault values..../5h init,a1,14h      iv(25) =,i3)
  230  if (iv(dtype) == iv(dtype0)) go to 250
           if (m == 0) write(pu,260) which
           m = 1
           write(pu,240) iv(dtype)
  240     format(20h dtype      iv(16) =,i3)
  250  i = 1
        j = jlim(alg)
        k = epslon
        l = iv(parsav)
        ndfalt = ndflt(alg)
      
        do ii = 1, ndfalt
           if (v(k) == v(l)) go to 280
                if (m == 0) write(pu,260) which
  260          format(/1h ,3a4,9halues..../)
                m = 1
                write(pu,270) vn(1,i), vn(2,i), k, v(k)
  270          format(1x,2a4,5h.. v(,i2,3h) =,d15.7)
  280     k = k + 1
           l = l + 1
           i = i + 1
           if (i == j) i = ijmp
        end do
      
        iv(dtype0) = iv(dtype)
        parsv1 = iv(parsav)
        call vcopy(iv(nvdflt), v(parsv1), v(epslon))
        return
      
  300  iv(1) = 15
        if (pu == 0) then
          return
        end if
        write(pu,310) liv, miv2
       310  format(/10h /// liv =,i5,17h must be at least,i5)
        if (liv < miv1) then
          return
        end if
        if (lv < iv(lastv)) go to 320
        return
      
  320  iv(1) = 16
        if (pu == 0) then
          return
        end if
        write(pu,330) lv, iv(lastv)
  330  format(/9h /// lv =,i5,17h must be at least,i5)
        return
      
  340  iv(1) = 67
        if (pu == 0) then
          return
        end if
        write(pu,350) alg
  350  format(/10h /// alg =,i5,15h must be 1 or 2)
      
        return
      end
      function reldst ( p, d, x, x0 )
      
      !*****************************************************************************80
      !
      !! RELDST computes the relative difference between X and X0.
      !
        integer p
      
        real ( kind = 8 ) reldst
        real ( kind = 8 ) d(p), x(p), x0(p)
        integer i
        real ( kind = 8 ) emax, t, xmax
      
        emax = 0.0D+00
        xmax = 0.0D+00
      
        do i = 1, p
          t = abs(d(i) * (x(i) - x0(i)))
          if (emax < t) emax = t
          t = d(i) * (abs(x(i)) + abs(x0(i)))
          if (xmax < t) xmax = t
        end do
      
        reldst = 0.0D+00
        if ( xmax > 0.0D+00 ) reldst = emax / xmax
      
        return
      end
      function rmdcon ( k )
      
      !*****************************************************************************80
      !
      !! RMDCON returns machine dependent constants.
      !
      !  Discussion:
      !
      !    Comments below contain data statements for various machines.
      !    To convert to another machine, place a c in column 1 of the
      !    data statement line(s) that correspond to the current machine
      !    and remove the c from column 1 of the data statement line(s)
      !    that correspond to the new machine.
      !
      !    the constant returned depends on k...
      !
      !         k = 1... smallest pos. eta such that -eta exists.
      !         k = 2... square root of eta.
      !         k = 3... unit roundoff = smallest pos. no. machep such
      !                  that 1 + machep > 1 .and. 1 - machep < 1.
      !         k = 4... square root of machep.
      !         k = 5... square root of big (see k = 6).
      !         k = 6... largest machine no. big such that -big exists.
      !
        integer k
        real ( kind = 8 ) rmdcon
        real ( kind = 8 ) big, eta, machep
        integer bigi(4), etai(4), machei(4)
        equivalence (big,bigi(1)), (eta,etai(1)), (machep,machei(1))
      !
      !  ibm 360, ibm 370, or xerox
      !
      !     data big/z7fffffffffffffff/, eta/z0010000000000000/,
      !    1     machep/z3410000000000000/
      !
      !  data general
      !
      !     data big/0.7237005577d+76/, eta/0.5397605347d-78/,
      !    1     machep/2.22044605d-16/
      !
      !  dec 11
      !
      !     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/
      !
      !  hp3000
      !
      !     data big/1.157920892d+77/, eta/8.636168556d-78/,
      !    1     machep/5.551115124d-17/
      !
      !  honeywell
      !
      !     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/
      !
      !  dec10
      !
      !     data big/"377777100000000000000000/,
      !    1     eta/"002400400000000000000000/,
      !    2     machep/"104400000000000000000000/
      !
      !  burroughs
      !
      !     data big/o0777777777777777,o7777777777777777/,
      !    1     eta/o1771000000000000,o7770000000000000/,
      !    2     machep/o1451000000000000,o0000000000000000/
      !
      !  control data
      !
      !     data big/37767777777777777777b,37167777777777777777b/,
      !    1     eta/00014000000000000000b,00000000000000000000b/,
      !    2     machep/15614000000000000000b,15010000000000000000b/
      !
      !  prime
      !
      !     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/
      !
      !  univac
      !
      !     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/
      !
      !  vax
      !
        data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/
      !
      !  cray 1
      !
      !     data bigi(1)/577767777777777777777b/,
      !    1     bigi(2)/000007777777777777776b/,
      !    2     etai(1)/200004000000000000000b/,
      !    3     etai(2)/000000000000000000000b/,
      !    4     machei(1)/377224000000000000000b/,
      !    5     machei(2)/000000000000000000000b/
      !
      !  port library -- requires more than just a data statement...
      !
      !     external d1mach
      !     real ( kind = 8 ) d1mach, zero
      !     data big/0.d+0/, eta/0.d+0/, machep/0.d+0/, zero/0.d+0/
      !     if (big > 0.0D+00) go to 1
      !        big = d1mach(2)
      !        eta = d1mach(1)
      !        machep = d1mach(4)
      !1    continue
      !
      ! end of port
      !
      !  body -
      !
        go to (10, 20, 30, 40, 50, 60), k
      
   10   rmdcon = eta
        return
      
   20   rmdcon = sqrt(256.d+0*eta)/16.d+0
        return
      
   30   rmdcon = machep
        return
      
   40   rmdcon = sqrt(machep)
        return
      
   50   rmdcon = sqrt(big/256.d+0)*16.d+0
        return
      
   60   rmdcon = big
      
        return
      end
      subroutine sgrad2 ( alpha, d, eta0, fx, g, irc, n, w, x )
      
      !*****************************************************************************80
      !
      !! SGRAD2 computes finite difference gradient by Stewart's scheme.
      !
      !  Discussion:
      !
      !    This subroutine uses an embellished form of the finite difference
      !    scheme proposed by Stewart to approximate the gradient of the
      !    function f(x), whose values are supplied by reverse communication.
      !
      !  Reference:
      !
      !    G W Stewart,
      !    A Modification of Davidon's Minimization Method to Accept
      !    Difference
      !    Approximations of Derivatives,
      !    Journal of the Association for Computing Machinery,
      !    Volume 14, pages. 72-83, 1967.
      !
      !  Parameters:
      !
      !  alpha in  (approximate) diagonal elements of the hessian of f(x).
      !      d in  scale vector such that d(i)*x(i), i = 1,...,n, are in
      !             comparable units.
      !   eta0 in  estimated bound on relative error in the function value...
      !             (true value) = (computed value)*(1+e),   where
      !             abs(e) <= eta0.
      !     fx i/o on input,  fx  must be the computed value of f(x).  on
      !             output with irc = 0, fx has been restored to its original
      !             value, the one it had when sgrad2 was last called with
      !             irc = 0.
      !      g i/o on input with irc = 0, g should contain an approximation
      !             to the gradient of f near x, e.g., the gradient at the
      !             previous iterate.  when sgrad2 returns with irc = 0, g is
      !             the desired finite-difference approximation to the
      !             gradient at x.
      !    irc i/o input/return code... before the very first call on sgrad2,
      !             the caller must set irc to 0.  whenever sgrad2 returns a
      !             nonzero value for irc, it has perturbed some component of
      !             x... the caller should evaluate f(x) and call sgrad2
      !             again with fx = f(x).
      !      n in  the number of variables (components of x) on which f
      !             depends.
      !      x i/o on input with irc = 0, x is the point at which the
      !             gradient of f is desired.  on output with irc nonzero, x
      !             is the point at which f should be evaluated.  on output
      !             with irc = 0, x has been restored to its original value
      !             (the one it had when sgrad2 was last called with irc = 0)
      !             and g contains the desired gradient approximation.
      !      w i/o work vector of length 6 in which sgrad2 saves certain
      !             quantities while the caller is evaluating f(x) at a
      !             perturbed x.
      !
      !      application and usage restrictions
      !
      !        this routine is intended for use with quasi-newton routines
      !     for unconstrained minimization (in which case  alpha  comes from
      !     the diagonal of the quasi-newton hessian approximation).
      !
      !      algorithm notes
      !
      !        this code departs from the scheme proposed by stewart (ref. 1)
      !     in its guarding against overly large or small step sizes and its
      !     handling of special cases (such as zero components of alpha or g).
      !
        integer irc, n
        real ( kind = 8 ) alpha(n), d(n), eta0, fx, g(n), w(6), x(n)
        external rmdcon
        real ( kind = 8 ) rmdcon
        integer fh, fx0, hsave, i, xisave
        real ( kind = 8 ) aai, afx, afxeta, agi, alphai, axi, axibar
        real ( kind = 8 ) discon, eta, gi, h, hmin
        real ( kind = 8 ) c2000, four, hmax0, hmin0, h0, machep, one, &
                         p002
        real ( kind = 8 ) three, two
      
        parameter (c2000=2.0d+3, four=4.0d+0, hmax0=0.02d+0, hmin0=5.0d+1 )
        parameter ( one=1.0d+0, p002=0.002d+0, three=3.0d+0 )
        parameter ( two=2.0d+0 )
      
        parameter (fh=3, fx0=4, hsave=5, xisave=6)
      !
        if (irc) 140, 100, 210
      !
      !      fresh start -- get machine-dependent constants
      !
      !     store machep in w(1) and h0 in w(2), where machep is the unit
      !     roundoff (the smallest positive number such that
      !     1 + machep > 1  and  1 - machep < 1),  and  h0 is the
      !     square-root of machep.
      !
  100  continue
      
        w(1) = rmdcon(3)
        w(2) = sqrt(w(1))
        w(fx0) = fx
      !
      !  increment  i  and start computing  g(i)
      !
  110  continue
      
           i = iabs(irc) + 1
        if (i > n) go to 300
           irc = i
           afx = abs(w(fx0))
           machep = w(1)
           h0 = w(2)
           hmin = hmin0 * machep
           w(xisave) = x(i)
           axi = abs(x(i))
           axibar = max (axi, one/d(i))
           gi = g(i)
           agi = abs(gi)
           eta = abs(eta0)
           if (afx > 0.0D+00) eta = max (eta, agi*axi*machep/afx)
           alphai = alpha(i)
           if (alphai == 0.0D+00) go to 170
           if (gi == 0.0D+00 .or. fx == 0.0D+00) go to 180
           afxeta = afx*eta
           aai = abs(alphai)
      !
      !  compute h = stewart's forward-difference step size.
      !
           if (gi**2 <= afxeta*aai) go to 120
                h = two*sqrt(afxeta/aai)
                h = h*(one - aai*h/(three*aai*h + four*agi))
                go to 130
      
  120 continue
      
           h = two*(afxeta*agi/(aai**2))**(one/three)
           h = h*(one - two*agi/(three*aai*h + four*agi))
      !
      !  ensure that  h  is not insignificantly small
      !
  130 continue
      
           h = max (h, hmin*axibar)
      !
      !  use forward difference if bound on truncation error is at
      !  most 10**-3.
      !
           if (aai*h <= p002*agi) go to 160
      !
      !  compute h = stewart*s step for central difference.
      !
           discon = c2000*afxeta
           h = discon/(agi + sqrt(gi**2 + aai*discon))
      !
      !  Ensure that  h  is neither too small nor too big
      !
           h = max (h, hmin*axibar)
           if (h >= hmax0*axibar) h = axibar * h0**(two/three)
      !
      !  Compute central difference
      !
           irc = -i
           go to 200
      
  140 continue
      
           h = -w(hsave)
           i = iabs(irc)
           if (h > 0.0D+00) go to 150
           w(fh) = fx
           go to 200
      
  150 continue
      
           g(i) = (w(fh) - fx) / (two * h)
           x(i) = w(xisave)
           go to 110
      !
      !  Compute forward differences in various cases
      !
  160 continue
      
           if (h >= hmax0*axibar) h = h0 * axibar
           if (alphai*gi < 0.0D+00) h = -h
           go to 200
      
  170 continue
           h = axibar
           go to 200
      
  180 continue
      
           h = h0 * axibar
      
  200 continue
      
           x(i) = w(xisave) + h
           w(hsave) = h
           return
      !
      !  compute actual forward difference
      !
  210 continue
      
           g(irc) = (fx - w(fx0)) / w(hsave)
           x(irc) = w(xisave)
           go to 110
      !
      !  Restore fx and indicate that g has been computed
      !
  300  continue
      
        fx = w(fx0)
        irc = 0
      
        return
      end
      subroutine slvmul ( p, y, s, x )
      
      !*****************************************************************************80
      !
      !! SLVMUL sets y = S * x.
      !
      !  Discussion:
      !
      !    s = p x p symmetric matrix.
      !    lower triangle of  s  stored rowwise.
      !
        integer p
        real ( kind = 8 ) s(*), x(p), y(p)
      !     dimension s(p*(p+1)/2)
        integer i, im1, j, k
        real ( kind = 8 ) xi
        real ( kind = 8 ) dotprd
      
        j = 1
      
        do i = 1, p
          y(i) = dotprd(i, s(j), x)
          j = j + i
        end do
      
        if (p <= 1) then
          return
        end if
      
        j = 1
      
        do i = 2, p
      
          xi = x(i)
          im1 = i - 1
          j = j + 1
      
          do k = 1, im1
            y(k) = y(k) + s(j)*xi
            j = j + 1
          end do
      
        end do
      
        return
      end
      subroutine smsno ( n, d, x, calcf, iv, liv, lv, v, uiparm, &
      urparm,ufparm )
      
      !*****************************************************************************80
      !
      !! SMSNO minimizes a general unconstrained objective function.
      !
      !  Discussion:
      !
      !    The routine uses finite-difference gradients and secant hessian
      !    approximations.
      !
      !    This routine interacts with SNOIT in an attempt
      !    to find an n-vector  x*  that minimizes the (unconstrained)
      !    objective function computed by  calcf.  (often the  x*  found is
      !    a local minimizer rather than a global one.)
      !
      !  Reference:
      !
      !    G W Stewart,
      !    A Modification of Davidon's Minimization Method to Accept
      !    Difference
      !      Approximations of Derivatives,
      !    Journal of the Association for Computing Machinery,
      !    Volume 14, pages 72-83, 1967.
      !
      !  Parameters:
      !
      !        the parameters for smsno are the same as those for sumsl
      !     (which see), except that calcg is omitted.  instead of calling
      !     calcg to obtain the gradient of the objective function at x,
      !     smsno calls sgrad2, which computes an approximation to the
      !     gradient by finite (forward and central) differences using the
      !     method of ref. 1.  the following input component is of interest
      !     in this regard (and is not described in sumsl).
      !
      ! v(eta0)  v(42) is an estimated bound on the relative error in the
      !             objective function value computed by calcf...
      !                  (true value) = (computed value) * (1 + e),
      !             where abs(e) <= v(eta0).  default = machep * 10**3,
      !             where machep is the unit roundoff.
      !
      !        the output values iv(nfcall) and iv(ngcall) have different
      !     meanings for smsno than for sumsl...
      !
      ! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
      !             function evaluations) excluding those made only for
      !             computing gradients.  the input value iv(mxfcal) is a
      !             limit on iv(nfcall).
      ! iv(ngcall)... iv(30) is the number of function evaluations made only
      !             for computing gradients.  the total number of function
      !             evaluations is thus  iv(nfcall) + iv(ngcall).
      !

        !---- Lines added by Andy ----
        use m_optimization
        use m_meamparameters

        logical firstEval
        !-----------------------------

        integer n, liv, lv
        integer iv(liv), uiparm(*)
        real ( kind = 8 ) d(n), x(n), v(lv), urparm(*)
      !     dimension v(77 + n*(n+17)/2), uiparm(*), urparm(*)
        external calcf, ufparm
        integer nf
        real ( kind = 8 ) fx
        integer nfcall, toobig
        parameter (nfcall=6, toobig=2)
     
        !---- Lines added by Andy ----
        firstEval=.true.
        opt_failed=.false.
        !-----------------------------
        do
      
          !print *,'calling snoit'
          call snoit ( d, fx, iv, liv, lv, n, v, x )
          !print *,'finished call...'

          if ( iv(1) > 2 ) then
            exit
          end if
      
          nf = iv(nfcall)
      
          !print *,'calling calcf'
          call calcf ( n, x, nf, fx, uiparm, urparm, ufparm )
          !print *,'finished call...'
          !---- Lines Added by Andy ----
          if (firstEval) then
             if (fx.lt.lowestoptfuncrand) then
                lowestoptfuncrand=fx
             endif
             if (n_optfunc.le.noptfuncstore) then
               if (fx.gt.upplimoptfunc) then
                  opt_failed=.true.
                  write(*,'(A58,G20.12,A17,G20.12,A2,G20.12,A1)') &
        ' Random potential parameters sampling... (current optfunc:', &
                      fx,', lowest optfunc:',lowestoptfuncrand,' >', &
                      upplimoptfunc,')'
                  return
               else
                  write(6,'(A57,G20.12,A2,G20.12,A40)') &
        ' Random potential parameters sampling... (lowest optfunc:', &
                       fx,' <',upplimoptfunc, &
               ')                                       '
               endif
             else
               if (fx.gt.upplimoptfunc_GA) then
                  opt_failed=.true.
                  write(*,'(A58,G20.12,A17,G20.12,A2,G20.12,A1)') &
        ' Random potential parameters sampling... (current optfunc:', &
                      fx,', lowest optfunc:',lowestoptfuncrand,' >', &
                      upplimoptfunc_GA,')'
                  return
               else
                  write(6,'(A57,G20.12,A2,G20.12,A40)') &
        ' Random potential parameters sampling... (lowest optfunc:', &
                       fx,' <',upplimoptfunc_GA, &
               ')                                       '
               endif
             endif
             firstEval=.false.
          !else
          !  !Check if optfunc is less than ... after ... fn evals.
          !  if ((nf.gt.n_optdiffint).and.(v(10).gt.optdiffint)) then
          !    !Andy write
          ! write(*,'(A9,F15.10,A14,F15.10,A7,I5,A26)') &
          !       ' Optfunc=',fx,' greater than ',optdiffint, &
          !       ' after ',nf,' function evals: Stopping.'
          !     opt_failed=.true.
          !     return
          !  endif
          endif
         !if (opt_failed.eqv..true.) then
         !   !If v(f) was found larger than optdiffint after nf
         !   !function evals (in the 'itsum' subroutine) then also
         !   !exit the CG optimizer.
         !   return
         !endif
          !-----------------------------

 
          if ( nf <= 0 ) then
            iv(toobig) = 1
          end if
      
        end do
      
        return
      end
      subroutine snoit ( d, fx, iv, liv, lv, n, v, x )
      
      !*****************************************************************************80
      !
      !! SNOIT is the iteration driver for SMSNO.
      !
      !  Discussion:
      !
      !    This routine minimizes a general unconstrained objective function
      !    using
      !    finite-difference gradients and secant hessian approximations.
      !
      !    This routine interacts with subroutine  sumit  in an attempt
      !    to find an n-vector  x*  that minimizes the (unconstrained)
      !    objective function  fx = f(x)  computed by the caller.  (often
      !    the  x*  found is a local minimizer rather than a global one.)
      !
      !  Reference:
      !
      !    G W Stewart,
      !    A Modification of Davidon's Minimization Method to Accept
      !    Difference
      !    Approximations of Derivatives,
      !    Journal of the Association for Computing Machinery,
      !    Volume 14, pages. 72-83, 1967.
      !
      !  Parameters:
      !
      !        the parameters for snoit are the same as those for sumsl
      !     (which see), except that calcf, calcg, uiparm, urparm, and ufparm
      !     are omitted, and a parameter  fx  for the objective function
      !     value at x is added.  instead of calling calcg to obtain the
      !     gradient of the objective function at x, snoit calls sgrad2,
      !     which computes an approximation to the gradient by finite
      !     (forward and central) differences using the method of ref. 1.
      !     the following input component is of interest in this regard
      !     (and is not described in sumsl).
      !
      ! v(eta0)  v(42) is an estimated bound on the relative error in the
      !             objective function value computed by calcf...
      !                  (true value) = (computed value) * (1 + e),
      !             where abs(e) <= v(eta0).  default = machep * 10**3,
      !             where machep is the unit roundoff.
      !
      !        the output values iv(nfcall) and iv(ngcall) have different
      !     meanings for smsno than for sumsl...
      !
      ! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
      !             function evaluations) excluding those made only for
      !             computing gradients.  the input value iv(mxfcal) is a
      !             limit on iv(nfcall).
      ! iv(ngcall)... iv(30) is the number of function evaluations made only
      !             for computing gradients.  the total number of function
      !             evaluations is thus  iv(nfcall) + iv(ngcall).
      !
        integer liv, lv, n
        integer iv(liv)
        real ( kind = 8 ) d(n), fx, x(n), v(lv)
      !     dimension v(77 + n*(n+17)/2)
      !
        external deflt, dotprd, sgrad2, sumit, vscopy
        real ( kind = 8 ) dotprd
        integer alpha, g1, i, iv1, j, k, w
        real ( kind = 8 ) zero
      
        integer eta0, f, g, lmat, nextv, nfgcal, ngcall
        integer niter, sgirc, toobig, vneed
      
        parameter ( eta0=42, f=10, g=28, lmat=42, nextv=47 )
        parameter ( nfgcal=7, ngcall=30, niter=31, sgirc=57 )
        parameter ( toobig=2, vneed=4)
      
        parameter ( zero=0.d+0)
      
        iv1 = iv(1)
        if (iv1 == 1) go to 10
        if (iv1 == 2) go to 50
        if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
        iv1 = iv(1)
        if (iv1 == 12 .or. iv1 == 13) iv(vneed) = iv(vneed) + 2*n + 6
        if (iv1 == 14) go to 10
        if (iv1 > 2 .and. iv1 < 12) go to 10
        g1 = 1
        if (iv1 == 12) iv(1) = 13
        go to 20
      
  10   g1 = iv(g)
      
  20   call sumit(d, fx, v(g1), iv, liv, lv, n, v, x)
        if (iv(1) - 2) 999, 30, 70
      !
      !  Compute gradient
      !
  30   if (iv(niter) == 0) call vscopy(n, v(g1), zero)
        j = iv(lmat)
        k = g1 - n
      
        do i = 1, n
          v(k) = dotprd(i, v(j), v(j))
          k = k + 1
          j = j + i
        end do
      !
      !  Undo increment of iv(ngcall) done by sumit
      !
        iv(ngcall) = iv(ngcall) - 1
      !
      !  Store return code from sgrad2 in iv(sgirc)
      !
        iv(sgirc) = 0
      !
      !  x may have been restored, so copy back fx...
      !
        fx = v(f)
        go to 60
      !
      !  gradient loop
      !
  50   if (iv(toobig) == 0) go to 60
        iv(nfgcal) = 0
        go to 10
      
  60   g1 = iv(g)
        alpha = g1 - n
        w = alpha - 6
        call sgrad2(v(alpha), d, v(eta0), fx, v(g1), iv(sgirc), n, &
             v(w),x)
        if (iv(sgirc) == 0) go to 10
           iv(ngcall) = iv(ngcall) + 1
           return
      
  70   if (iv(1) /= 14) then
              return
            end if
      !
      !  Storage allocation
      !
        iv(g) = iv(nextv) + n + 6
        iv(nextv) = iv(g) + n
        if (iv1 /= 13) go to 10
      
       999  continue
      
        return
      end
      function stopx ( )
      
      !*****************************************************************************80
      !
      !! STOPX checks to see if the BREAK key has been pressed.
      !
      !  Discussion:
      !
      !     this function may serve as the stopx (asynchronous interruption)
      !     function for the nl2sol (nonlinear least-squares) package at
      !     those installations which do not wish to implement a
      !     dynamic stopx.
      !
      !     at installations where the nl2sol system is used
      !     interactively, this dummy stopx should be replaced by a
      !     function that returns .true. if and only if the interrupt
      !     (break) key has been pressed since the last call on stopx.
      !
        logical stopx
      
        stopx = .false.
      
        return
      end
      subroutine sumit ( d, fx, g, iv, liv, lv, n, v, x)
      
      !*****************************************************************************80
      !
      !! SUMIT carries out unconstrained minimization iterations for SUMSL.
      !
      !  Discussion:
      !
      !    The routine uses double-dogleg/BFGS steps.
      !
      !    parameters iv, n, v, and x are the same as the corresponding
      !    ones to sumsl (which see), except that v can be shorter (since
      !    the part of v that sumsl uses for storing g is not needed).
      !    moreover, compared with sumsl, iv(1) may have the two additional
      !    output values 1 and 2, which are explained below, as is the use
      !    of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
      !    output value from sumsl (and smsno), is not referenced by
      !    sumit or the subroutines it calls.
      !
      !    fx and g need not have been initialized when sumit is called
      !    with iv(1) = 12, 13, or 14.
      !
      ! iv(1) = 1 means the caller should set fx to f(x), the function value
      !             at x, and call sumit again, having changed none of the
      !             other parameters.  an exception occurs if f(x) cannot be
      !             (e.g. if overflow would occur), which may happen because
      !             of an oversized step.  in this case the caller should set
      !             iv(toobig) = iv(2) to 1, which will cause sumit to ig-
      !             nore fx and try a smaller step.  the parameter nf that
      !             sumsl passes to calcf (for possible use by calcg) is a
      !             copy of iv(nfcall) = iv(6).
      ! iv(1) = 2 means the caller should set g to g(x), the gradient vector
      !             of f at x, and call sumit again, having changed none of
      !             the other parameters except possibly the scale vector d
      !             when iv(dtype) = 0.  the parameter nf that sumsl passes
      !             to calcg is iv(nfgcal) = iv(7).  if g(x) cannot be
      !             evaluated, then the caller may set iv(nfgcal) to 0, in
      !             which case sumit will return with iv(1) = 65.
      !
      !  Parameters:
      !
      ! d.... scale vector.
      ! fx... function value.
      ! g.... gradient vector.
      ! iv... integer value array.
      ! liv.. length of iv (at least 60).
      ! lv... length of v (at least 71 + n*(n+13)/2).
      ! n.... number of variables (components in x and g).
      ! v.... floating-point value array.
      ! x.... vector of parameters to be optimized.
      !
        integer liv
        integer lv
        integer n
      
        integer iv(liv)
        real ( kind = 8 ) d(n)
        real ( kind = 8 ) fx
        real ( kind = 8 ) g(n)
        real ( kind = 8 ) v(lv)
        real ( kind = 8 ) x(n)
        integer dg1, g01, i, k, l, lstgst, nwtst1, step1
        integer        temp1, w, x01, z
        real ( kind = 8 ) t
        real ( kind = 8 ) half, negone, one, onep2, zero
        logical stopx
        real ( kind = 8 ) dotprd, reldst, v2norm
        integer cnvcod, dg, dgnorm, dinit, dstnrm, dst0, f, f0, fdif
        integer gthg, gtstep, g0, incfac, inith, irc, kagqt, lmat, lmax0
        integer lmaxs, mode, model, mxfcal, mxiter, nextv, nfcall, &
                nfgcal
        integer ngcall, niter, nreduc, nwtstp, preduc, radfac, radinc
        integer radius, rad0, reldx, restor, step, stglim, stlstg, &
                toobig
        integer tuner4, tuner5, vneed, xirc, x0
      
        parameter (cnvcod=55, dg=37, g0=48, inith=25, irc=29, kagqt=33 )
        parameter ( mode=35, model=5, mxfcal=17, mxiter=18, nfcall=6 )
        parameter ( nfgcal=7, ngcall=30, niter=31, nwtstp=34, radinc=8 )
        parameter ( restor=9, step=40, stglim=11, stlstg=41, toobig=2 )
        parameter ( vneed=4, xirc=13, x0=43)
      
        parameter (dgnorm=1, dinit=38, dstnrm=2, dst0=3, f=10, f0=13 )
        parameter ( fdif=11, gthg=44, gtstep=4, incfac=23, lmat=42 )
        parameter ( lmax0=35, lmaxs=36, nextv=47, nreduc=6, preduc=7 )
        parameter ( radfac=16, radius=8, rad0=9, reldx=17, tuner4=29 )
        parameter ( tuner5=30)
      
        parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, &
         onep2=1.2d+0,zero=0.d+0)
      !
        i = iv(1)
        if (i == 1) go to 50
        if (i == 2) go to 60
      !
      !   check validity of iv and v input values
      !
        if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
        if (iv(1) == 12 .or. iv(1) == 13) then
          iv(vneed) = iv(vneed) + n*(n+13)/2
        end if
        call parck(2, d, iv, liv, lv, n, v)
        i = iv(1) - 2
        if (i > 12) then
          return
        end if
        go to (180, 180, 180, 180, 180, 180, 120, 90, 120, 10, 10, 20),&
           i
      !
      !   storage allocation
      !
   10    l = iv(lmat)
        iv(x0) = l + n*(n+1)/2
        iv(step) = iv(x0) + n
        iv(stlstg) = iv(step) + n
        iv(g0) = iv(stlstg) + n
        iv(nwtstp) = iv(g0) + n
        iv(dg) = iv(nwtstp) + n
        iv(nextv) = iv(dg) + n
        if (iv(1) /= 13) go to 20
           iv(1) = 14
           return
      !
      !   initialization
      !
   20   iv(niter) = 0
        iv(nfcall) = 1
        iv(ngcall) = 1
        iv(nfgcal) = 1
        iv(mode) = -1
        iv(model) = 1
        iv(stglim) = 1
        iv(toobig) = 0
        iv(cnvcod) = 0
        iv(radinc) = 0
        v(rad0) = 0.0D+00
        if (v(dinit) >= 0.0D+00) call vscopy(n, d, v(dinit))
        if (iv(inith) /= 1) go to 40
      !
      !  set the initial hessian approximation to diag(d)**-2
      !
           l = iv(lmat)
           call vscopy(n*(n+1)/2, v(l), zero)
           k = l - 1
      
           do i = 1, n
             k = k + i
             t = d(i)
             if (t <= 0.0D+00) t = one
             v(k) = t
           end do
      !
      !  compute initial function value
      !
   40   iv(1) = 1
        return
      
   50   v(f) = fx
        if (iv(mode) >= 0) go to 180
        iv(1) = 2
        if (iv(toobig) == 0) then
          return
        end if
           iv(1) = 63
           go to 300
      !
      !   make sure gradient could be computed
      !
   60   if (iv(nfgcal) /= 0) go to 70
           iv(1) = 65
           go to 300
      
   70   dg1 = iv(dg)
        call vvmulp(n, v(dg1), g, d, -1)
        v(dgnorm) = v2norm(n, v(dg1))
      
        if (iv(cnvcod) /= 0) go to 290
        if (iv(mode) == 0) go to 250
      !
      !   allow first step to have scaled 2-norm at most v(lmax0)
      !
        v(radius) = v(lmax0)
      
        iv(mode) = 0
      !
      !  main loop
      !
      !   print iteration summary, check iteration limit
      !
   80   call itsum(d, g, iv, liv, lv, n, v, x)
   90   k = iv(niter)
        if (k < iv(mxiter)) go to 100
           iv(1) = 10
           go to 300
      !
      !   update radius
      !
  100  iv(niter) = k + 1
        if(k>0)v(radius) = v(radfac) * v(dstnrm)
      !
      !   initialize for start of next iteration
      !
        g01 = iv(g0)
        x01 = iv(x0)
        v(f0) = v(f)
        iv(irc) = 4
        iv(kagqt) = -1
      !
      !      copy x to x0, g to g0
      !
        call vcopy(n, v(x01), x)
        call vcopy(n, v(g01), g)
      !
      !  Check STOPX and function evaluation limit
      !
  110  if ( .not. stopx ( ) ) go to 130
           iv(1) = 11
           go to 140
      !
      !  Come here when restarting after func. eval. limit or STOPX.
      !
  120  if (v(f) >= v(f0)) go to 130
           v(radfac) = one
           k = iv(niter)
           go to 100
      
  130  if (iv(nfcall) < iv(mxfcal)) go to 150
           iv(1) = 9
  140     if (v(f) >= v(f0)) go to 300
      !
      !  in case of STOPX or function evaluation limit with
      !  improved v(f), evaluate the gradient at x.
      !
                iv(cnvcod) = iv(1)
                go to 240
      !
      !  Compute candidate step
      !
  150  step1 = iv(step)
        dg1 = iv(dg)
        nwtst1 = iv(nwtstp)
        if (iv(kagqt) >= 0) go to 160
           l = iv(lmat)
           call livmul(n, v(nwtst1), v(l), g)
           v(nreduc) = half * dotprd(n, v(nwtst1), v(nwtst1))
           call litvmu(n, v(nwtst1), v(l), v(nwtst1))
           call vvmulp(n, v(step1), v(nwtst1), d, 1)
           v(dst0) = v2norm(n, v(step1))
           call vvmulp(n, v(dg1), v(dg1), d, -1)
           call ltvmul(n, v(step1), v(l), v(dg1))
           v(gthg) = v2norm(n, v(step1))
           iv(kagqt) = 0
  160  call dbdog(v(dg1), lv, n, v(nwtst1), v(step1), v)
        if (iv(irc) == 6) go to 180
      !
      !   check whether evaluating f(x0 + step) looks worthwhile
      !
        if (v(dstnrm) <= 0.0D+00) go to 180
        if (iv(irc) /= 5) go to 170
        if (v(radfac) <= one) go to 170
        if (v(preduc) <= onep2 * v(fdif)) go to 180
      !
      !  Compute f(x0 + step)
      !
  170  x01 = iv(x0)
        step1 = iv(step)
        call vaxpy(n, x, one, v(step1), v(x01))
        iv(nfcall) = iv(nfcall) + 1
        iv(1) = 1
        iv(toobig) = 0
        return
      !
      !  Assess candidate step.
      !
  180  x01 = iv(x0)
        v(reldx) = reldst(n, d, x, v(x01))
        call assst(iv, liv, lv, v)
        step1 = iv(step)
        lstgst = iv(stlstg)
        if (iv(restor) == 1) call vcopy(n, x, v(x01))
        if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
        if (iv(restor) /= 3) go to 190
           call vcopy(n, v(step1), v(lstgst))
           call vaxpy(n, x, one, v(step1), v(x01))
           v(reldx) = reldst(n, d, x, v(x01))
      
  190  k = iv(irc)
        go to (200,230,230,230,200,210,220,220,220,220,220,220,280, &
               250), k
      !
      !      recompute step with changed radius
      !
  200     v(radius) = v(radfac) * v(dstnrm)
           go to 110
      !
      !   compute step of length v(lmaxs) for singular convergence test.
      !
  210  v(radius) = v(lmaxs)
        go to 150
      !
      !   convergence or false convergence
      !
  220  iv(cnvcod) = k - 4
        if (v(f) >= v(f0)) go to 290
           if (iv(xirc) == 14) go to 290
                iv(xirc) = 14
      !
      !  Process acceptable step.
      !
  230  if (iv(irc) /= 3) go to 240
           step1 = iv(step)
           temp1 = iv(stlstg)
      !
      !      set  temp1 = hessian * step  for use in gradient tests
      !
           l = iv(lmat)
           call ltvmul(n, v(temp1), v(l), v(step1))
           call lvmul(n, v(temp1), v(l), v(temp1))
      !
      !   compute gradient
      !
  240  iv(ngcall) = iv(ngcall) + 1
        iv(1) = 2
        return
      !
      !   initializations -- g0 = g - g0, etc.
      !
  250  g01 = iv(g0)
        call vaxpy(n, v(g01), negone, v(g01), g)
        step1 = iv(step)
        temp1 = iv(stlstg)
        if (iv(irc) /= 3) go to 270
      !
      !   set v(radfac) by gradient tests
      !
      !  Set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))
      !
           call vaxpy(n, v(temp1), negone, v(g01), v(temp1))
           call vvmulp(n, v(temp1), v(temp1), d, -1)
      !
      !  Do gradient tests
      !
           if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) then
             go to 260
           end if
      
           if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  then
             go to 270
           end if
      
  260               v(radfac) = v(incfac)
      !
      !   update h, loop
      !
  270  w = iv(nwtstp)
        z = iv(x0)
        l = iv(lmat)
        call wzbfgs(v(l), n, v(step1), v(w), v(g01), v(z))
      !
      !  Use the n-vectors starting at v(step1) and v(g01) for scratch.
      !
        call lupdat(v(temp1), v(step1), v(l), v(g01), v(l), n, v(w), &
        v(z))
        iv(1) = 2
        go to 80
      !
      !   misc. details
      !
      !   bad parameters to assess
      !
  280  iv(1) = 64
        go to 300
      !
      !  Print summary of final iteration and other requested items
      !
  290  iv(1) = iv(cnvcod)
        iv(cnvcod) = 0
  300  call itsum(d, g, iv, liv, lv, n, v, x)
      
        return
      end
      subroutine sumsl(n, d, x, calcf, calcg, iv, liv, lv, v, uiparm, &
          urparm,ufparm)
      
      !*****************************************************************************80
      !
      !! SUMSL minimizes a general unconstrained objective function.
      !
      !  Discussion:
      !
      !    The routine uses analytic gradient and hessian approximation from
      !    the secant update.
      !
      !    This routine interacts with subroutine  sumit  in an attempt
      !    to find an n-vector  x*  that minimizes the (unconstrained)
      !    objective function computed by  calcf.  (often the  x*  found is
      !    a local minimizer rather than a global one.)
      !
      !  Reference:
      !
      !    J E Dennis, David Gay, and R E Welsch,
      !    An Adaptive Nonlinear Least-squares Algorithm,
      !    ACM Transactions on Mathematical Software,
      !    Volume 7, Number 3, 1981.
      !
      !    J E Dennis, H H W Mei,
      !    Two New Unconstrained Optimization Algorithms Which Use
      !    Function and Gradient Values,
      !    Journal of Optimization Theory and Applications,
      !    Volume 28, pages 453-482, 1979.
      !
      !    J E Dennis, Jorge More,
      !    Quasi-Newton Methods, Motivation and Theory,
      !    SIAM Review,
      !    Volume 19, pages 46-89, 1977.
      !
      !    D Goldfarb,
      !    Factorized Variable Metric Methods for Unconstrained Optimization,
      !    Mathematics of Computation,
      !    Volume 30, pages 796-811, 1976.
      !
      !  Parameters:
      !
      ! n  (input) the number of variables on which  f  depends, i.e.,
      !                  the number of components in  x.
      ! d  (input/output) a scale vector such that  d(i)*x(i),
      !                  i = 1,2,...,n,  are all in comparable units.
      !                  d can strongly affect the behavior of sumsl.
      !                  finding the best choice of d is generally a trial-
      !                  and-error process.  choosing d so that d(i)*x(i)
      !                  has about the same value for all i often works well.
      !                  the defaults provided by subroutine deflt (see iv
      !                  below) require the caller to supply d.
      ! x........ (input/output) before (initially) calling sumsl, the call-
      !                  er should set  x  to an initial guess at  x*.  when
      !                  sumsl returns,  x  contains the best point so far
      !                  found, i.e., the one that gives the least value so
      !                  far seen for  f(x).
      ! calcf.... (input) a subroutine that, given x, computes f(x).  calcf
      !                  must be declared external in the calling program.
      !                  it is invoked by
      !                       call calcf(n, x, nf, f, uiparm, urparm, ufparm)
      !                  when calcf is called, nf is the invocation
      !                  count for calcf.  nf is included for possible use
      !                  with calcg.  if x is out of bounds (e.g., if it
      !                  would cause overflow in computing f(x)), then calcf
      !                  should set nf to 0.  this will cause a shorter step
      !                  to be attempted.  (if x is in bounds, then calcf
      !                  should not change nf.)  the other parameters are as
      !                  described above and below.  calcf should not change
      !                  n, p, or x.
      ! calcg.... (input) a subroutine that, given x, computes g(x), the gra-
      !                  dient of f at x.  calcg must be declared external in
      !                  the calling program.  it is invoked by
      !                       call calcg(n, x, nf, g, uiparm, urparm, ufaprm)
      !                  when calcg is called, nf is the invocation
      !                  count for calcf at the time f(x) was evaluated.  the
      !                  x passed to calcg is usually the one passed to calcf
      !                  on either its most recent invocation or the one
      !                  prior to it.  if calcf saves intermediate results
      !                  for use by calcg, then it is possible to tell from
      !                  nf whether they are valid for the current x (or
      !                  which copy is valid if two copies are kept).  if g
      !                  cannot be computed at x, then calcg should set nf to
      !                  0.  in this case, sumsl will return with iv(1) = 65.
      !                  (if g can be computed at x, then calcg should not
      !                  changed nf.)  the other parameters to calcg are as
      !                  described above and below.  calcg should not change
      !                  n or x.
      ! iv....... (input/output) an integer value array of length liv (see
      !                  below) that helps control the sumsl algorithm and
      !                  that is used to store various intermediate quanti-
      !                  ties.  of particular interest are the initialization/
      !                  return code iv(1) and the entries in iv that control
      !                  printing and limit the number of iterations and func-
      !                  tion evaluations.  see the section on iv input
      !                  values below.
      ! liv...... (input) length of iv array.  must be at least 60.  if liv
      !                  is too small, then sumsl returns with iv(1) = 15.
      !                  when sumsl returns, the smallest allowed value of
      !                  liv is stored in iv(lastiv) -- see the section on
      !                  iv output values below.  (this is intended for use
      !                  with extensions of sumsl that handle constraints.)
      ! lv....... (input) length of v array.  must be at least 71+n*(n+15)/2.
      !                  (at least 77+n*(n+17)/2 for smsno, at least
      !                  78+n*(n+12) for humsl).  if lv is too small, then
      !                  sumsl returns with iv(1) = 16.  when sumsl returns,
      !                  the smallest allowed value of lv is stored in
      !                  iv(lastv) -- see the section on iv output values
      !                  below.
      ! v........ (input/output) a floating-point value array of length lv
      !                  (see below) that helps control the sumsl algorithm
      !                  and that is used to store various intermediate
      !                  quantities.  of particular interest are the entries
      !                  in v that limit the length of the first step
      !                  attempted (lmax0) and specify convergence tolerances
      !                  (afctol, lmaxs, rfctol, sctol, xctol, xftol).
      ! uiparm... (input) user integer parameter array passed without change
      !                  to calcf and calcg.
      ! urparm... (input) user floating-point parameter array passed without
      !                  change to calcf and calcg.
      ! ufparm... (input) user external subroutine or function passed without
      !                  change to calcf and calcg.
      !
      !   iv input values (from subroutine deflt)
      !
      ! iv(1)...  on input, iv(1) should have a value between 0 and 14......
      !             0 and 12 mean this is a fresh start.  0 means that
      !                  deflt(2, iv, liv, lv, v)
      !             is to be called to provide all default values to iv and
      !             v.  12 (the value that deflt assigns to iv(1)) means the
      !             caller has already called deflt and has possibly changed
      !             some iv and/or v entries to non-default values.
      !             13 means deflt has been called and that sumsl (and
      !             sumit) should only do their storage allocation.  that is,
      !             they should set the output components of iv that tell
      !             where various subarrays arrays of v begin, such as iv(g)
      !             (and, for humsl and humit only, iv(dtol)), and return.
      !             14 means that a storage has been allocated (by a call
      !             with iv(1) = 13) and that the algorithm should be
      !             started.  when called with iv(1) = 13, sumsl returns
      !             iv(1) = 14 unless liv or lv is too small (or n is not
      !             positive).  default = 12.
      ! iv(inith).... iv(25) tells whether the hessian approximation h should
      !             be initialized.  1 (the default) means sumit should
      !             initialize h to the diagonal matrix whose i-th diagonal
      !             element is d(i)**2.  0 means the caller has supplied a
      !             cholesky factor  l  of the initial hessian approximation
      !             h = l*(l**t)  in v, starting at v(iv(lmat)) = v(iv(42))
      !             (and stored compactly by rows).  note that iv(lmat) may
      !             be initialized by calling sumsl with iv(1) = 13 (see
      !             the iv(1) discussion above).  default = 1.
      ! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
      !             (calls on calcf) allowed.  if this number does not suf-
      !             fice, then sumsl returns with iv(1) = 9.  default = 200.
      ! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
      !             it also indirectly limits the number of gradient evalua-
      !             tions (calls on calcg) to iv(mxiter) + 1.  if iv(mxiter)
      !             iterations do not suffice, then sumsl returns with
      !             iv(1) = 10.  default = 150.
      ! iv(outlev)... iv(19) controls the number and length of iteration sum-
      !             mary lines printed (by itsum).  iv(outlev) = 0 means do
      !             not print any summary lines.  otherwise, print a summary
      !             line after each abs(iv(outlev)) iterations.  if iv(outlev)
      !             is positive, then summary lines of length 78 (plus carri-
      !             age control) are printed, including the following...  the
      !             iteration and function evaluation counts, f = the current
      !             function value, relative difference in function values
      !             achieved by the latest step (i.e., reldf = (f0-v(f))/f01,
      !             where f01 is the maximum of abs(v(f)) and abs(v(f0)) and
      !             v(f0) is the function value from the previous itera-
      !             tion), the relative function reduction predicted for the
      !             step just taken (i.e., preldf = v(preduc) / f01, where
      !             v(preduc) is described below), the scaled relative change
      !             in x (see v(reldx) below), the step parameter for the
      !             step just taken (stppar = 0 means a full newton step,
      !             between 0 and 1 means a relaxed newton step, between 1
      !             and 2 means a double dogleg step, greater than 2 means
      !             a scaled down Cauchy step -- see subroutine dbldog), the
      !             2-norm of the scale vector d times the step just taken
      !             (see v(dstnrm) below), and npreldf, i.e.,
      !             v(nreduc)/f01, where v(nreduc) is described below -- if
      !             npreldf is positive, then it is the relative function
      !             reduction predicted for a newton step (one with
      !             stppar = 0).  if npreldf is negative, then it is the
      !             negative of the relative function reduction predicted
      !             for a step computed with step bound v(lmaxs) for use in
      !             testing for singular convergence.
      !                  if iv(outlev) is negative, then lines of length 50
      !             are printed, including only the first 6 items listed
      !             above (through reldx).
      !             default = 1.
      ! iv(parprt)... iv(20) = 1 means print any nondefault v values on a
      !             fresh start or any changed v values on a restart.
      !             iv(parprt) = 0 means skip this printing.  default = 1.
      ! iv(prunit)... iv(21) is the output unit number on which all printing
      !             is done.  iv(prunit) = 0 means suppress all printing.
      !             default = standard output unit (unit 6 on most systems).
      ! iv(solprt)... iv(22) = 1 means print out the value of x returned (as
      !             well as the gradient and the scale vector d).
      !             iv(solprt) = 0 means skip this printing.  default = 1.
      ! iv(statpr)... iv(23) = 1 means print summary statistics upon return-
      !             ing.  these consist of the function value, the scaled
      !             relative change in x caused by the most recent step (see
      !             v(reldx) below), the number of function and gradient
      !             evaluations (calls on calcf and calcg), and the relative
      !             function reductions predicted for the last step taken and
      !             for a newton step (or perhaps a step bounded by v(lmaxs)
      !             -- see the descriptions of preldf and npreldf under
      !             iv(outlev) above).
      !             iv(statpr) = 0 means skip this printing.
      !             iv(statpr) = -1 means skip this printing as well as that
      !             of the one-line termination reason message.  default = 1.
      ! iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
      !             (on a fresh start only).  iv(x0prt) = 0 means skip this
      !             printing.  default = 1.
      !
      !   (selected) iv output values
      !
      ! iv(1)........ on output, iv(1) is a return code....
      !             3 = x-convergence.  the scaled relative difference (see
      !                  v(reldx)) between the current parameter vector x and
      !                  a locally optimal parameter vector is very likely at
      !                  most v(xctol).
      !             4 = relative function convergence.  the relative differ-
      !                  ence between the current function value and its lo-
      !                  cally optimal value is very likely at most v(rfctol).
      !             5 = both x- and relative function convergence (i.e., the
      !                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
      !             6 = absolute function convergence.  the current function
      !                  value is at most v(afctol) in absolute value.
      !             7 = singular convergence.  the hessian near the current
      !                  iterate appears to be singular or nearly so, and a
      !                  step of length at most v(lmaxs) is unlikely to yield
      !                  a relative function decrease of more than v(sctol).
      !             8 = false convergence.  the iterates appear to be converg-
      !                  ing to a noncritical point.  this may mean that the
      !                  convergence tolerances (v(afctol), v(rfctol),
      !                  v(xctol)) are too small for the accuracy to which
      !                  the function and gradient are being computed, that
      !                  there is an error in computing the gradient, or that
      !                  the function or gradient is discontinuous near x.
      !             9 = function evaluation limit reached without other con-
      !                  vergence (see iv(mxfcal)).
      !            10 = iteration limit reached without other convergence
      !                  (see iv(mxiter)).
      !            11 = STOPX returned .true. (external interrupt).  see the
      !                  usage notes below.
      !            14 = storage has been allocated (after a call with
      !                  iv(1) = 13).
      !            17 = restart attempted with n changed.
      !            18 = d has a negative component and iv(dtype) <= 0.
      !            19...43 = v(iv(1)) is out of range.
      !            63 = f(x) cannot be computed at the initial x.
      !            64 = bad parameters passed to assess (which should not
      !                  occur).
      !            65 = the gradient could not be computed at x (see calcg
      !                  above).
      !            67 = bad first parameter to deflt.
      !            80 = iv(1) was out of range.
      !            81 = n is not positive.
      ! iv(g)........ iv(28) is the starting subscript in v of the current
      !             gradient vector (the one corresponding to x).
      ! iv(lastiv)... iv(44) is the least acceptable value of liv.  (it is
      !             only set if liv is at least 44.)
      ! iv(lastv).... iv(45) is the least acceptable value of lv.  (it is
      !             only set if liv is large enough, at least iv(lastiv).)
      ! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
      !             function evaluations).
      ! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
      !             calcg).
      ! iv(niter).... iv(31) is the number of iterations performed.
      !
      !   (selected) v input values (from subroutine deflt)
      !
      ! v(bias)..... v(43) is the bias parameter used in subroutine dbldog --
      !             see that subroutine for details.  default = 0.8.
      ! v(afctol)... v(31) is the absolute function convergence tolerance.
      !             if sumsl finds a point where the function value is less
      !             than v(afctol) in absolute value, and if sumsl does not
      !             return with iv(1) = 3, 4, or 5, then it returns with
      !             iv(1) = 6.  this test can be turned off by setting
      !             v(afctol) to zero.  default = max(10**-20, machep**2),
      !             where machep is the unit roundoff.
      ! v(dinit).... v(38), if nonnegative, is the value to which the scale
      !             vector d is initialized.  default = -1.
      ! v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
      !             very first step that sumsl attempts.  this parameter can
      !             markedly affect the performance of sumsl.
      ! v(lmaxs).... v(36) is used in testing for singular convergence -- if
      !             the function reduction predicted for a step of length
      !             bounded by v(lmaxs) is at most v(sctol) * abs(f0), where
      !             f0  is the function value at the start of the current
      !             iteration, and if sumsl does not return with iv(1) = 3,
      !             4, 5, or 6, then it returns with iv(1) = 7.  default = 1.
      ! v(rfctol)... v(32) is the relative function convergence tolerance.
      !             if the current model predicts a maximum possible function
      !             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0)
      !             at the start of the current iteration, where  f0  is the
      !             then current function value, and if the last step attempt-
      !             ed achieved no more than twice the predicted function
      !             decrease, then sumsl returns with iv(1) = 4 (or 5).
      !             default = max(10**-10, machep**(2/3)), where machep is
      !             the unit roundoff.
      ! v(sctol).... v(37) is the singular convergence tolerance -- see the
      !             description of v(lmaxs) above.
      ! v(tuner1)... v(26) helps decide when to check for false convergence.
      !             this is done if the actual function decrease from the
      !             current step is no more than v(tuner1) times its predict-
      !             ed value.  default = 0.1.
      ! v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
      !             (see v(nreduc)) is tried that has v(reldx) <= v(xctol)
      !             and if this step yields at most twice the predicted func-
      !             tion decrease, then sumsl returns with iv(1) = 3 (or 5).
      !             (see the description of v(reldx) below.)
      !             default = machep**0.5, where machep is the unit roundoff.
      ! v(xftol).... v(34) is the false convergence tolerance.  if a step is
      !             tried that gives no more than v(tuner1) times the predict-
      !             ed function decrease and that has v(reldx) <= v(xftol),
      !             and if sumsl does not return with iv(1) = 3, 4, 5, 6, or
      !             7, then it returns with iv(1) = 8.  (see the description
      !             of v(reldx) below.)  default = 100*machep, where
      !             machep is the unit roundoff.
      ! v(*)........ deflt supplies to v a number of tuning constants, with
      !             which it should ordinarily be unnecessary to tinker.  see
      !             section 17 of version 2.2 of the nl2sol usage summary
      !             (i.e., the appendix to ref. 1) for details on v(i),
      !             i = decfac, incfac, phmnfc, phmxfc, rdfcmn, rdfcmx,
      !             tuner2, tuner3, tuner4, tuner5.
      !
      !   (selected) v output values
      !
      ! v(dgnorm)... v(1) is the 2-norm of (diag(d)**-1)*g, where g is the
      !             most recently computed gradient.
      ! v(dstnrm)... v(2) is the 2-norm of diag(d)*step, where step is the
      !             current step.
      ! v(f)........ v(10) is the current function value.
      ! v(f0)....... v(13) is the function value at the start of the current
      !             iteration.
      ! v(nreduc)... v(6), if positive, is the maximum function reduction
      !             possible according to the current model, i.e., the func-
      !             tion reduction predicted for a newton step (i.e.,
      !             step = -h**-1 * g,  where  g  is the current gradient and
      !             h is the current hessian approximation).
      !                  if v(nreduc) is negative, then it is the negative of
      !             the function reduction predicted for a step computed with
      !             a step bound of v(lmaxs) for use in testing for singular
      !             convergence.
      ! v(preduc)... v(7) is the function reduction predicted (by the current
      !             quadratic model) for the current step.  this (divided by
      !             v(f0)) is used in testing for relative function
      !             convergence.
      ! v(reldx).... v(17) is the scaled relative change in x caused by the
      !             current step, computed as
      !                  max(abs(d(i)*(x(i)-x0(i)), 1 <= i <= p) /
      !                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p),
      !             where x = x0 + step.
      !
      !  notes
      !
      !   algorithm notes
      !
      !        this routine uses a hessian approximation computed from the
      !     bfgs update (see ref 3).  only a cholesky factor of the hessian
      !     approximation is stored, and this is updated using ideas from
      !     ref. 4.  steps are computed by the double dogleg scheme described
      !     in ref. 2.  the steps are assessed as in ref. 1.
      !
      !   usage notes
      !
      !        after a return with iv(1) <= 11, it is possible to restart,
      !     i.e., to change some of the iv and v input values described above
      !     and continue the algorithm from the point where it was interrupt-
      !     ed.  iv(1) should not be changed, nor should any entries of iv
      !     and v other than the input values (those supplied by deflt).
      !        those who do not wish to write a calcg which computes the
      !     gradient analytically should call smsno rather than sumsl.
      !     smsno uses finite differences to compute an approximate gradient.
      !        those who would prefer to provide f and g (the function and
      !     gradient) by reverse communication rather than by writing subrou-
      !     tines calcf and calcg may call on sumit directly.  see the com-
      !     ments at the beginning of sumit.
      !        those who use sumsl interactively may wish to supply their
      !     own STOPX function, which should return .true. if the break key
      !     has been pressed since STOPX was last invoked.  this makes it
      !     possible to externally interrupt sumsl (which will return with
      !     iv(1) = 11 if STOPX returns .true.).
      !        storage for g is allocated at the end of v.  thus the caller
      !     may make v longer than specified above and may allow calcg to use
      !     elements of g beyond the first n as scratch storage.
      !
      !   portability notes
      !
      !        the sumsl distribution tape contains both single- and double-
      !     precision versions of the sumsl source code, so it should be un-
      !     necessary to change precisions.
      !        only the function rmdcon contains machine-dependent
      !     constants.  to change from one machine to another, it should
      !     suffice to change the (few) relevant lines in these functions.
      !        intrinsic functions are explicitly declared.  on certain com-
      !     puters (e.g. univac), it may be necessary to comment out these
      !     declarations.  so that this may be done automatically by a simple
      !     program, such declarations are preceded by a comment having c/+
      !     in columns 1-3 and blanks in columns 4-72 and are followed by
      !     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
      !        the sumsl source code is expressed in 1966 ansi standard
      !     fortran.  it may be converted to fortran 77 by commenting out all
      !     lines that fall between a line having c/6 in columns 1-3 and a
      !     line having c/7 in columns 1-3 and by removing (i.e., replacing
      !     by a blank) the c in column 1 of the lines that follow the c/7
      !     line and precede a line having c/ in columns 1-2 and blanks in
      !     columns 3-72.  these changes convert some data statements into
      !     parameter statements, convert some variables from real to
      !     character*4, and make the data statements that initialize these
      !     variables use character strings delimited by primes instead
      !     of hollerith constants.  (such variables and data statements
      !     appear only in modules itsum and parck.  parameter statements
      !     appear nearly everywhere.)  these changes also add save state-
      !     ments for variables given machine-dependent constants by rmdcon.
      !
        integer n, liv, lv
        integer iv(liv), uiparm(*)
        real ( kind = 8 ) d(n), x(n), v(lv), urparm(*)
      !     dimension v(71 + n*(n+15)/2), uiparm(*), urparm(*)
      
        integer g1, iv1, nf
        real ( kind = 8 ) f
        integer nextv, nfcall, nfgcal, g, toobig, vneed
      
        parameter (nextv=47, nfcall=6, nfgcal=7, g=28, toobig=2, &
            vneed=4)
      
        external ufparm
      
        if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
        iv1 = iv(1)
        if (iv1 == 12 .or. iv1 == 13) iv(vneed) = iv(vneed) + n
        if (iv1 == 14) go to 10
        if (iv1 > 2 .and. iv1 < 12) go to 10
        g1 = 1
        if (iv1 == 12) iv(1) = 13
        go to 20
      
   10   g1 = iv(g)
      
   20   call sumit(d, f, v(g1), iv, liv, lv, n, v, x)
        if (iv(1) - 2) 30, 40, 50
      
   30   nf = iv(nfcall)
        call calcf(n, x, nf, f, uiparm, urparm, ufparm)
        if (nf <= 0) iv(toobig) = 1
        go to 20
      
   40   call calcg(n, x, iv(nfgcal), v(g1), uiparm, urparm, ufparm)
        go to 20
      
   50   if (iv(1) /= 14) then
              return
            end if
      !
      !  Storage allocation
      !
        iv(g) = iv(nextv)
        iv(nextv) = iv(g) + n
        if (iv1 /= 13) go to 10
      
        return
      end
      subroutine timestamp ( )
      
      !*****************************************************************************80
      !
      !! TIMESTAMP prints the current YMDHMS date as a time stamp.
      !
      !  Example:
      !
      !    31 May 2001   9:45:54.872 AM
      !
      !  Modified:
      !
      !    06 August 2005
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    None
      !
        implicit none
      
        character ( len = 8 ) ampm
        integer ( kind = 4 ) d
        integer ( kind = 4 ) h
        integer ( kind = 4 ) m
        integer ( kind = 4 ) mm
        character ( len = 9 ), parameter, dimension(12) :: month = (/ &
          'January  ', 'February ', 'March    ', 'April    ', &
          'May      ', 'June     ', 'July     ', 'August   ', &
          'September', 'October  ', 'November ', 'December ' /)
        integer ( kind = 4 ) n
        integer ( kind = 4 ) s
        integer ( kind = 4 ) values(8)
        integer ( kind = 4 ) y
      
        call date_and_time ( values = values )
      
        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)
      
        if ( h < 12 ) then
          ampm = 'AM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
          else
            ampm = 'PM'
          end if
        else
          h = h - 12
          if ( h < 12 ) then
            ampm = 'PM'
          else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
              ampm = 'Midnight'
            else
              ampm = 'AM'
            end if
          end if
        end if
      
        write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3, &
                1x,a)' ) d, trim ( month(m) ), y, h, ':', n, ':', &
                s, '.', mm, trim ( ampm )
      
        return
      end
      function v2norm ( p, x )
      
      !*****************************************************************************80
      !
      !! V2NORM returns the 2-norm of the p-vector X.
      !
      !  Discussion:
      !
      !    The routine tries to avoid underflow.
      !
      !  Parameters:
      !
        integer p
      
        real ( kind = 8 ) x(p)
        integer i, j
        real ( kind = 8 ) r, scale
        real ( kind = 8 ), save :: sqteta = 0.0D+00
        real ( kind = 8 ) t, xi
        real ( kind = 8 ) rmdcon
        real ( kind = 8 ) v2norm
      
        v2norm = 0.0D+00
      
        if (p <= 0 ) then
          return
        end if
      
        if ( all ( x(1:p) == 0.0D+00 ) ) then
          return
        end if
      
        scale = 0.0D+00
        do i = 1, p
          if ( x(i) /= 0.0D+00 ) then
            scale = abs(x(i))
            exit
          end if
        end do
      
        if ( scale == 0.0D+00 ) then
          return
        end if
      
        if ( p <= i ) then
          v2norm = scale
          return
        end if
      
        t = 1.0D+00
        if ( sqteta == 0.0D+00 ) then
          sqteta = rmdcon(2)
        end if
      !
      !  sqteta is (slightly larger than) the square root of the
      !  smallest positive floating point number on the machine.
      !  the tests involving sqteta are done to prevent underflows.
      !
        j = i + 1
        do i = j, p
          xi = abs(x(i))
          if (xi <= scale) then
            r = xi / scale
            if (r > sqteta) t = t + r*r
          else
            r = scale / xi
            if (r <= sqteta) r = 0.0D+00
            t = 1.0D+00  +  t * r*r
            scale = xi
          end if
        end do
      
        v2norm = scale * sqrt(t)
      
        return
      end
      subroutine vaxpy ( p, w, a, x, y )
      
      !*****************************************************************************80
      !
      !! VAXPY sets w = a*x + y.
      !
      !  Discussion:
      !
      !    w, x, y = p-vectors, a = scalar
      !
      !  Parameters:
      !
        implicit none
      
        integer p
      
        real ( kind = 8 ) a
        real ( kind = 8 ) w(p)
        real ( kind = 8 ) x(p)
        real ( kind = 8 ) y(p)
      
        w(1:p) = a * x(1:p) + y(1:p)
      
        return
      end
      subroutine vcopy ( p, y, x )
      
      !*****************************************************************************80
      !
      !! VCOPY sets y = x.
      !
      !  Discussion:
      !
      !    x and y are p-vectors
      !
        implicit none
      
        integer p
        real ( kind = 8 ) x(p)
        real ( kind = 8 ) y(p)
      
        y(1:p) = x(1:p)
      
        return
      end
      subroutine vdflt ( alg, lv, v )
      
      !*****************************************************************************80
      !
      !! VDFLT supplies default values to V.
      !
      !  Discussion:
      !
      !    alg = 1 means regression constants.
      !    alg = 2 means general unconstrained optimization constants.
      !
        implicit none
      
        integer alg, lv
        real ( kind = 8 ) v(lv)
        real ( kind = 8 ) rmdcon
        real ( kind = 8 ) machep, mepcrt, one, sqteps, three
        integer afctol, bias, cosmin, decfac, delta0, dfac, dinit, &
                dltfdc
        integer dltfdj, dtinit, d0init, epslon, eta0, fuzz, huberc
        integer incfac, lmax0, lmaxs, phmnfc, phmxfc, rdfcmn, rdfcmx
        integer rfctol, rlimit, rsptol, sctol, sigmin, tuner1, tuner2
        integer tuner3, tuner4, tuner5, xctol, xftol
      
        parameter (one=1.d+0, three=3.d+0)
      
        parameter (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44 )
        parameter ( dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39 )
        parameter ( d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48 )
        parameter ( incfac=23, lmax0=35, lmaxs=36, phmnfc=20, &
                    phmxfc=21 )
        parameter ( rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49 )
        parameter ( sctol=37, sigmin=50, tuner1=26, tuner2=27, &
                    tuner3=28 )
        parameter ( tuner4=29, tuner5=30, xctol=33, xftol=34)
      
        machep = rmdcon(3)
        v(afctol) = 1.d-20
      
        if ( machep > 1.d-10 ) then
          v(afctol) = machep**2
        end if
      
        v(decfac) = 0.5d+0
        sqteps = rmdcon(4)
        v(dfac) = 0.6d+0
        v(delta0) = sqteps
        v(dtinit) = 1.d-6
        mepcrt = machep ** (one/three)
        v(d0init) = 1.d+0
        v(epslon) = 0.1d+0
        v(incfac) = 2.d+0
        v(lmax0) = 1.d+0
        v(lmaxs) = 1.d+0
        v(phmnfc) = -0.1d+0
        v(phmxfc) = 0.1d+0
        v(rdfcmn) = 0.1d+0
        v(rdfcmx) = 4.d+0
        v(rfctol) = max (1.d-10, mepcrt**2)
        v(sctol) = v(rfctol)
        v(tuner1) = 0.1d+0
        v(tuner2) = 1.d-4
        v(tuner3) = 0.75d+0
        v(tuner4) = 0.5d+0
        v(tuner5) = 0.75d+0
        v(xctol) = sqteps
        v(xftol) = 1.d+2 * machep
      
        if ( alg < 2 ) then
          v(cosmin) = max (1.d-6, 1.d+2 * machep)
          v(dinit) = 0.d+0
          v(dltfdc) = mepcrt
          v(dltfdj) = sqteps
          v(fuzz) = 1.5d+0
          v(huberc) = 0.7d+0
          v(rlimit) = rmdcon(5)
          v(rsptol) = 1.d-3
          v(sigmin) = 1.d-4
        else
          v(bias) = 0.8d+0
          v(dinit) = -1.0d+0
          v(eta0) = 1.0d+3 * machep
        end if
      
        return
      end
      subroutine vscopy ( p, y, s )
      
      !*****************************************************************************80
      !
      !! VSCOPY sets the vector Y to scalar S.
      !
        implicit none
      
        integer p
      
        real ( kind = 8 ) s
        real ( kind = 8 ) y(p)
      
        y(1:p) = s
      
        return
      end
      subroutine vvmulp ( n, x, y, z, k )
      
      !*****************************************************************************80
      !
      !! VVMULP sets x(i) = y(i) * z(i)**k, 1 <= i <= n (for k = 1 or -1)
      !
        implicit none
      
        integer n
      
        integer k
        real ( kind = 8 ) x(n)
        real ( kind = 8 ) y(n)
        real ( kind = 8 ) z(n)
      
        if ( k < 0 ) then
          x(1:n) = y(1:n) / z(1:n)
        else
          x(1:n) = y(1:n) * z(1:n)
        end if
      
        return
      end
      subroutine wzbfgs ( l, n, s, w, y, z )
      
      !*****************************************************************************80
      !
      !! WZBFGS compute Y and Z for LUPDAT corresponding to BFGS update.
      !
      !  Discussion:
      !
      !    When S is computed in certain ways, for example by GQTSTP or
      !    DBLDOG, it is possible to save N**2/2 operations since L'*S
      !    or L*L'*S is then known.
      !
      !    If the BFGS update to L*L' would reduce its determinant to
      !    less than EPS times its old value, then this routine in effect
      !    replaces Y by THETA*Y + (1-THETA)*L*L'*S, where THETA
      !    (between 0 and 1) is chosen to make the reduction factor = EPS.
      !
      !  Parameters:
      !
      !    l (i/o) cholesky factor of hessian, a lower triang. matrix stored
      !             compactly by rows.
      !
      !    n (input) order of  l  and length of  s,  w,  y,  z.
      !
      !    s (input) the step just taken.
      !
      !    w (output) right singular vector of rank 1 correction to l.
      !
      !    y (input) change in gradients corresponding to s.
      !
      !    z (output) left singular vector of rank 1 correction to l.
      !
        implicit none
      
        integer n
      
        real ( kind = 8 ) dotprd
        real ( kind = 8 ) cs
        real ( kind = 8 ) cy
        real ( kind = 8 ), parameter :: eps = 0.1D+00
        real ( kind = 8 ) epsrt
        real ( kind = 8 ) l(n*(n+1)/2)
        real ( kind = 8 ) s(n)
        real ( kind = 8 ) shs
        real ( kind = 8 ) theta
        real ( kind = 8 ) w(n)
        real ( kind = 8 ) y(n)
        real ( kind = 8 ) ys
        real ( kind = 8 ) z(n)
      
        call ltvmul ( n, w, l, s )
        shs = dotprd ( n, w, w )
        ys = dotprd ( n, y, s )
      
        if ( ys < eps * shs ) then
          theta = ( 1.0D+00 - eps ) * shs / ( shs - ys )
          epsrt = sqrt ( eps )
          cy = theta / ( shs * epsrt )
          cs = ( 1.0D+00 + ( theta - 1.0D+00 ) / epsrt ) / shs
        else
          cy = 1.0D+00 / ( sqrt ( ys ) * sqrt ( shs ) )
          cs = 1.0D+00 / shs
        end if
      
        call livmul ( n, z, l, y )
      
        z(1:n) = cy * z(1:n) - cs * w(1:n)
      
        return
      end
      
subroutine variables_to_p

    !---------------------------------------------------------------c
    !
    !     Copies the sensibly-named MEAM variables (cmin, etc)
    !     onto the p() array.
    !
    !     Called by:     program MEAMfit,initializemeam,meamenergy
    !     Calls:         -
    !     Arguments:     -
    !     Returns:       cmin,cmax,meamtau,
    !                 meamrhodecay,meame0,meamrho0,
    !                 meamemb3,meamemb4,pairpotparameter,
    !                 rs,rc
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff
    !
    !---------------------------------------------------------------c

    use m_filenames
    use m_meamparameters
    use m_optimization
    use m_atomproperties
    use m_generalinfo

    implicit none

    integer i,j,k,l,m

    j=1
    k=1
    l=1
    do i=1,m3
        p(i)=cmin(j,k,l)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        elseif (l.lt.maxspecies) then
            j=1
            k=1
            l=l+1
        endif
    enddo

    j=1
    k=1
    l=1
    do i=m3+1,2*m3
        p(i)=cmax(j,k,l)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        elseif (l.lt.maxspecies) then
            j=1
            k=1
            l=l+1
        endif
    enddo
    !      print *,'cmax=',cmax
    !      cmax(1:maxspecies,1:maxspecies,1:maxspecies)=p(1+m3:2*m3)

    j=1
    k=0
    do i=2*m3+1,2*m3+lm1*m1
        p(i)=meamtau(k,j)
        !print *,'p(',i,')=',p(i),', meamtau(',k,',',j,')=',meamtau(k,j)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.lmax) then
            j=1
            k=k+1
        endif
    enddo
    !print *,'in variables_to_p:'
    !print *
    !print *,'meamtau=',meamtau
    !print *,'p(2*m3+1,2*m3+lm1*m1)=',p(2*m3+1:2*m3+lm1*m1)
    !      meamtau(0:lmax,1:maxspecies)=p(2*m3+1:2*m3+lm1*m1)

    j=1
    k=0
    l=1
    m=1
    do i=2*m3+lm1*m1+1,2*m3+lm1*m1+12*lm1*m2
        p(i)=meamrhodecay(j,k,l,m)
        !print *,'conv to p(',i,')=',p(i)
        if (j.lt.12) then
            j=j+1
        elseif (k.lt.lmax) then
            j=1
            k=k+1
        elseif (m.lt.maxspecies) then
            j=1
            k=0
            m=m+1
        elseif (l.lt.maxspecies) then
            j=1
            k=0
            l=l+1
            m=1
        endif
    enddo
    p(2*m3+lm1*m1+12*lm1*m2+1: &
        2*m3+m1+lm1*m1+12*lm1*m2)=meame0(1:maxspecies)
    !      print *,'meame0=',meame0
    p(2*m3+m1+lm1*m1+12*lm1*m2+1: &
        2*m3+2*m1+lm1*m1+12*lm1*m2)=meamrho0(1:maxspecies)
    !      print *,'meamrho0=',meamrho0
    p(2*m3+2*m1+lm1*m1+12*lm1*m2+1: &
        2*m3+3*m1+lm1*m1+12*lm1*m2)=meamemb3(1:maxspecies)
    !      print *,'meamemb3=',meamemb3
    p(2*m3+3*m1+lm1*m1+12*lm1*m2+1: &
        2*m3+4*m1+lm1*m1+12*lm1*m2)=meamemb4(1:maxspecies)
    !      print *,'meamemb4=',meamemb4

    j=1
    if (maxspecies.eq.1) then
        do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,2*m3+(4+lm1)*m1+12*lm1*m2+32
            p(i)=pairpotparameter(j,1,1)
            j=j+1
        enddo
    else
        k=1
        l=1
        do i=2*m3+(4+lm1)*m1+12*lm1*m2+1,m4
            p(i)=pairpotparameter(j,k,l)
            if (j.lt.32) then
                j=j+1
            elseif (l.lt.maxspecies) then
                j=1
                l=l+1
            elseif (k.lt.maxspecies) then
                j=1
                l=1
                k=k+1
            endif
        enddo
    endif
    !      print *,'pairpotparameter=',pairpotparameter
    !      pairpotparameter(1:32,1,1)=p(2*m3+(4+lm1)*m1+6*lm1*m2+1:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+32)
    !      pairpotparameter(1:32,2,1)=p(2*m3+(4+lm1)*m1+6*lm1*m2+33:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+64)
    !      pairpotparameter(1:32,1,2)=p(2*m3+(4+lm1)*m1+6*lm1*m2+65:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+96)
    !      pairpotparameter(1:32,2,2)=p(2*m3+(4+lm1)*m1+6*lm1*m2+97:
    !     +             2*m3+(4+lm1)*m1+6*lm1*m2+128)

    j=1
    k=1
    do i=m4+1,m4+m2
        p(i)=rs(j,k)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        endif
    enddo

    j=1
    k=1
    do i=m4+m2+1,m4+2*m2
        p(i)=rc(j,k)
        if (j.lt.maxspecies) then
            j=j+1
        elseif (k.lt.maxspecies) then
            j=1
            k=k+1
        endif
    enddo

    j=1
    do i=m4+2*m2+1,m4+2*m2+m1
        p(i)=enconst(j)
        j=j+1
    enddo

end subroutine variables_to_p
subroutine writeBoundsForRandParas

    !--------------------------------------------------------------c
    !
    !     Reads in the limits which are to be used to initialize
    !     the random starting values of the MEAM parameters from
    !     the 'dataforMEAMparagen' file.
    !
    !     Called by:
    !     Calls:
    !     Arguments:
    !     Returns:
    !     Files read:
    !     Files written:
    !
    !     Andy Duff
    !
    !--------------------------------------------------------------c

    use m_meamparameters
    use m_filenames
    use m_atomproperties
    use m_generalinfo
    use m_optimization

    implicit none

    integer i,j,ii,ll,ispc,jspc,spc_cnt
    character*80 string
    !Variables just used for reading in start_meam_parameters file:
    integer nfreep

    open(unit=1,file=trim(settingsfile),access='append')

    write(1,*) "RandParasManual"

    !meamtau
    do ll=0,lmax
        do ispc=1,maxspecies
            if (freep(2*m3+ispc+ll*m1).ge.2) then
                write(1,*) meamtau_minorder(ll,ispc)
                write(1,*) meamtau_maxorder(ll,ispc)
            endif
        enddo
    enddo

    !electron density functions
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if ((ispc.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                do ll=0,lmax
                    !Check if we need to generate random paramaters for
                    !thi(ispc,jspc,ll)
                    nfreep=0
                    do i=1,12
                        if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+i).ge.2) then
                            nfreep=nfreep+1
                        endif
                    enddo
                    if (nfreep.gt.0) then
                        if (typethi(ispc,jspc).eq.1) then !Cubic form
                            !Read in parameters from 'dataforMEAMparagen'
                            !necessary
                            !to generate these parameters
                            write(1,*) meamrhodecay_maxradius(ll,ispc,jspc)
                            write(1,*) meamrhodecay_minradius(ll,ispc,jspc)
                            write(1,*) meamrhodecay_negvals(ll,ispc,jspc)

                            !For each coefficient which needs to be generated,
                            !read in
                            !min and max orders (10^..)
                            ii=1
                            do i=1,6
                                if (freep(2*m3+lm1*m1+12*lm1*spc_cnt+12*ll+2*i-1).ge.2) &
                                    then
                                    write(1,*) meamrhodecay_minorder(ii,ll,ispc,jspc)
                                    write(1,*) meamrhodecay_maxorder(ii,ll,ispc,jspc)
                                    ii=ii+1
                                endif
                            enddo
                        endif
                    endif
                enddo
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !embedding functions
    do ispc=1,maxspecies
        do i=1,4
            if (freep(2*m3+lm1*m1+12*lm1*m2+ispc+(i-1)*m1).ge.2) then
                if (i.eq.1) then
                    write(1,*) meame0_minorder(ispc)
                    write(1,*) meame0_maxorder(ispc)
                elseif (i.eq.2) then
                    write(1,*) meamrho0_minorder(ispc)
                    write(1,*) meamrho0_maxorder(ispc)
                elseif (i.eq.3) then
                    write(1,*) meamemb3_minorder(ispc)
                    write(1,*) meamemb3_maxorder(ispc)
                elseif (i.eq.4) then
                    write(1,*) meamemb4_minorder(ispc)
                    write(1,*) meamemb4_maxorder(ispc)
                endif
            endif
        enddo
    enddo

    !pair-potentials
    spc_cnt=0
    do ispc=1,maxspecies
        do jspc=1,maxspecies
            if (jspc.ge.ispc) then
                !Check if we need to generate random paramaters for
                !pairpot(isp,jspc)
                nfreeppairpot(ispc,jspc)=0
                do i=1,32
                    if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+i+32*spc_cnt) &
                        .ge.2) then
                        nfreeppairpot(ispc,jspc)=nfreeppairpot(ispc,jspc)+1
                    endif
                enddo
                if (nfreeppairpot(ispc,jspc).gt.0) then
                    if (typepairpot(ispc,jspc).eq.2) then
                        !Read in parameters from 'dataforMEAMparagen' necessary
                        !to
                        !generate these parameters
                        write(1,*) pairpotparameter_maxradius(ispc,jspc)
                        write(1,*) pairpotparameter_minradius(ispc,jspc)
                        write(1,*) pairpotparameter_negvals(ispc,jspc)
                        if ((pairpotparameter_negvals(ispc,jspc).lt.1).or. &
                            (pairpotparameter_negvals(ispc,jspc).gt.3)) then
                            print *,'      ERROR: you must specify 1, 2 or 3,'
                            print *,'      stopping.'
                            stop
                        endif
                        !For each coefficient which needs to be generated, read
                        !in min and max orders (10^..)
                        ii=1
                        do i=1,16
                            if (freep(2*m3+(4+lm1)*m1+12*lm1*m2+2*i-1+32*spc_cnt) &
                                .ge.2) then
                                write(1,*) pairpotparameter_minorder(ii,ispc,jspc)
                                write(1,*) pairpotparameter_maxorder(ii,ispc,jspc)
                                ii=ii+1
                            endif
                        enddo
                    else
                        print *,'typepairpot=',typepairpot(ispc,jspc), &
                            'not supported,stopping'
                        stop
                    endif
                endif
            endif
            spc_cnt=spc_cnt+1
        enddo
    enddo

    !enconst
    do ispc=1,maxspecies
        if (freep(m4+2*m2+ispc).ge.2) then
            write(1,*) enconst_minorder(ispc)
            write(1,*) enconst_maxorder(ispc)
        endif
    enddo

    close(1)

end subroutine writeBoundsForRandParas
subroutine writederivedquantities(outch)

    !---------------------------------------------------------------------c
    !
    !     Prints a set of user defined quantities derived from the
    !     energies (or potentially also the force components). This
    !     subroutine should be adjusted manually by the user to print out
    !     the desired quantities.
    !
    !     Called by:     program MEAMfit
    !     Calls:         -
    !     Arguments:     outch
    !     Returns:       -
    !     Files read:    -
    !     Files written: Depends on value of outch (6=standard output)
    !
    !     Andy Duff, Feb 2011.
    !
    !---------------------------------------------------------------------c

    use m_datapoints
    use m_optimization
    use m_geometry
    use m_filenames

    implicit none

    integer outch,i,idatapoint
    real(8) ECoct_fitted,ECtet_fitted,E2C100b_fitted,EFe_fitted, &
        ECoct_true,ECtet_true,E2C100b_true,EFe_true, &
        E2C111_fitted,E2C111_true,E2C1h0_fitted,E2C1h0_true, &
        E2C11h_fitted,E2C11h_true, &
        E2C110a_fitted,E2C110a_true,E2C110b_fitted,E2C110b_true, &
        EFeV_fitted,EFeV_true,E2CVcis_fitted,E2CVcis_true, &
        E2CVtrans_fitted,E2CVtrans_true, &
        E3CVfac_fitted,E3CVfac_true, &
        E3CVmer_fitted,E3CVmer_true, &
        ECV11h_fitted,ECV11h_true, &
        ECV1h0_fitted,ECV1h0_true, &
        ECV1hh_fitted,ECV1hh_true, &
        ECVh00_fitted,ECVh00_true, &
        ECVj00_fitted,ECVj00_true, &
        ECVhh0_fitted,ECVhh0_true, &
        ECVjh0_fitted,ECVjh0_true, &
        ECVjj0_fitted,ECVjj0_true, &
        ECVj1h_fitted,ECVj1h_true, &
        ECVsubs_fitted,ECVsubs_true, &
        ESI_fitted,ESI_true, &
        ECSIh00a_fitted,ECSIh00a_true

    ECoct_true=0d0
    ECtet_true=0d0
    E2C100b_true=0d0
    E2C111_true=0d0
    E2C1h0_true=0d0
    E2C11h_true=0d0
    E2C110a_true=0d0
    E2C110b_true=0d0
    ECVh00_true=0d0
    ECVj00_true=0d0
    ECVhh0_true=0d0
    ECV1h0_true=0d0
    ECV1hh_true=0d0
    ECV11h_true=0d0
    ECVjh0_true=0d0
    ECVjj0_true=0d0
    ECVj1h_true=0d0
    ECVsubs_true=0d0
    E2CVtrans_true=0d0
    E2CVcis_true=0d0
    E3CVfac_true=0d0
    E3CVmer_true=0d0
    EFe_true=0d0
    EFeV_true=0d0
    ESI_true=0d0
    ECSIh00a_true=0d0

    idatapoint=1
    do i=1,nstruct
        if (optimizeforce(i).ne.0) then
            idatapoint=idatapoint+3*gn_forces(i)
        elseif (optimizeforce(i).eq.0) then
            if (strucnames(i).eq.'Coct') then
                ECoct_fitted=fitdata(idatapoint)
                ECoct_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'Ctet') then
                ECtet_fitted=fitdata(idatapoint)
                ECtet_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2C100b') then
                E2C100b_fitted=fitdata(idatapoint)
                E2C100b_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2C1h0') then
                E2C1h0_fitted=fitdata(idatapoint)
                E2C1h0_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2C110a') then
                E2C110a_fitted=fitdata(idatapoint)
                E2C110a_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2C110b') then
                E2C110b_fitted=fitdata(idatapoint)
                E2C110b_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2C11h') then
                E2C11h_fitted=fitdata(idatapoint)
                E2C11h_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2C111') then
                E2C111_fitted=fitdata(idatapoint)
                E2C111_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'Fe') then
                EFe_fitted=fitdata(idatapoint)
                EFe_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'FeV') then
                EFeV_fitted=fitdata(idatapoint)
                EFeV_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2CVcis') then
                E2CVcis_fitted=fitdata(idatapoint)
                E2CVcis_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'2CVtrans') then
                E2CVtrans_fitted=fitdata(idatapoint)
                E2CVtrans_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'3CVfac') then
                E3CVfac_fitted=fitdata(idatapoint)
                E3CVfac_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'3CVmer') then
                E3CVmer_fitted=fitdata(idatapoint)
                E3CVmer_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CVh00') then
                ECVh00_fitted=fitdata(idatapoint)
                ECVh00_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CVj00') then
                ECVj00_fitted=fitdata(idatapoint)
                ECVj00_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CVhh0') then
                ECVhh0_fitted=fitdata(idatapoint)
                ECVhh0_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CVjj0') then
                ECVjj0_fitted=fitdata(idatapoint)
                ECVjj0_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CV11h') then
                ECV11h_fitted=fitdata(idatapoint)
                ECV11h_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CV1h0') then
                ECV1h0_fitted=fitdata(idatapoint)
                ECV1h0_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CV1hh') then
                ECV1hh_fitted=fitdata(idatapoint)
                ECV1hh_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CVjh0') then
                ECVjh0_fitted=fitdata(idatapoint)
                ECVjh0_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CVj1h') then
                ECVj1h_fitted=fitdata(idatapoint)
                ECVj1h_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CVsubs') then
                ECVsubs_fitted=fitdata(idatapoint)
                ECVsubs_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'SI') then
                ESI_fitted=fitdata(idatapoint)
                ESI_true=truedata(idatapoint)
            elseif (strucnames(i).eq.'CSIh00a') then
                ECSIh00a_fitted=fitdata(idatapoint)
                ECSIh00a_true=truedata(idatapoint)
            endif
            idatapoint=idatapoint+1
        endif
    enddo

    write(outch,*)
    write(outch,*) 'Derived quantities:'
    write(outch,*)
    write(outch,*) 'Structure          fitdata', &
        '       truedata'
    write(outch,*) '----------------------------------------', &
        '----'
    if ((ECoct_true.ne.0d0).and.(ECtet_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_migration       ', &
            ECtet_fitted-ECoct_fitted,'  ',ECtet_true-ECoct_true
    endif
    if ((ECoct_true.ne.0d0).and.(E2C100b_true.ne.0d0).and. &
        (EFe_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2C100b)       ', &
            E2C100b_fitted+EFe_fitted-2d0*ECoct_fitted,'  ', &
            E2C100b_true+EFe_true-2d0*ECoct_true
    endif
    if ((ECoct_true.ne.0d0).and.(E2C1h0_true.ne.0d0).and. &
        (EFe_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2C1h0)        ', &
            E2C1h0_fitted+EFe_fitted-2d0*ECoct_fitted,'  ', &
            E2C1h0_true+EFe_true-2d0*ECoct_true
    endif
    if ((ECoct_true.ne.0d0).and.(E2C110a_true.ne.0d0).and. &
        (EFe_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2C110a)       ', &
            E2C110a_fitted+EFe_fitted-2d0*ECoct_fitted,'  ', &
            E2C110a_true+EFe_true-2d0*ECoct_true
    endif
    if ((ECoct_true.ne.0d0).and.(E2C110b_true.ne.0d0).and. &
        (EFe_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2C110b)       ', &
            E2C110b_fitted+EFe_fitted-2d0*ECoct_fitted,'  ', &
            E2C110b_true+EFe_true-2d0*ECoct_true
    endif
    if ((ECoct_true.ne.0d0).and.(E2C11h_true.ne.0d0).and. &
        (EFe_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2C11h)        ', &
            E2C11h_fitted+EFe_fitted-2d0*ECoct_fitted,'  ', &
            E2C11h_true+EFe_true-2d0*ECoct_true
    endif
    if ((ECoct_true.ne.0d0).and.(E2C111_true.ne.0d0).and. &
        (EFe_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2C111)        ', &
            E2C111_fitted+EFe_fitted-2d0*ECoct_fitted,'  ', &
            E2C111_true+EFe_true-2d0*ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECVh00_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CVh00)        ', &
            ECVh00_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECVh00_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECVj00_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CVj00)        ', &
            ECVj00_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECVj00_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECVhh0_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CVhh0)        ', &
            ECVhh0_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECVhh0_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECVjj0_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CVjj0)        ', &
            ECVjj0_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECVjj0_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECV1h0_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CV1h0)        ', &
            ECV1h0_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECV1h0_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECV1hh_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CV1hh)        ', &
            ECV1hh_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECV1hh_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECV11h_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CV11h)        ', &
            ECV11h_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECV11h_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECVjh0_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CVjh0)        ', &
            ECVjh0_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECVjh0_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECVj1h_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CVj1h)        ', &
            ECVj1h_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECVj1h_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(EFe_true.ne.0d0) &
        .and.(ECVsubs_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CVsubs)       ', &
            ECVsubs_fitted+EFe_fitted-EFeV_fitted-ECoct_fitted,'  ', &
            ECVsubs_true+EFe_true-EFeV_true-ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(E2CVcis_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2CVcis)       ', &
            E2CVcis_fitted+2d0*EFe_fitted-EFeV_fitted-2d0*ECoct_fitted,'  ', &
            E2CVcis_true+2d0*EFe_true-EFeV_true-2d0*ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(E2CVtrans_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(2CVtrans)     ', &
            E2CVtrans_fitted+2d0*EFe_fitted-EFeV_fitted-2d0*ECoct_fitted,'  ' &
            ,E2CVtrans_true+2d0*EFe_true-EFeV_true-2d0*ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(E3CVmer_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(3CVmer)       ', &
            E3CVmer_fitted+3d0*EFe_fitted-EFeV_fitted-3d0*ECoct_fitted,'  ' &
            ,E3CVmer_true+3d0*EFe_true-EFeV_true-3d0*ECoct_true
    endif
    if ((EFeV_true.ne.0d0).and.(E3CVfac_true.ne.0d0).and. &
        (ECoct_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(3CVfac)       ', &
            E3CVfac_fitted+3d0*EFe_fitted-EFeV_fitted-3d0*ECoct_fitted,'  ' &
            ,E3CVfac_true+3d0*EFe_true-EFeV_true-3d0*ECoct_true
    endif
    if ((ESI_true.ne.0d0).and.(ECSIh00a_true.ne.0d0).and. &
        (ECoct_true.ne.0d0).and.(EFe_true.ne.0d0)) then
        write(outch,'(A19,F12.5,A2,F12.5)') 'E_b(CSIh00a)      ', &
            ECSIh00a_fitted+EFe_fitted-ESI_fitted-ECoct_fitted,'  ' &
            ,ECSIh00a_true+EFe_true-ESI_true-ECoct_true
    endif

end subroutine writederivedquantities
subroutine writemeamdatapoints(outch,extraOutput)

    !---------------------------------------------------------------------c
    !
    !     Prints meam and vasp energies and forces to the screen and to
    !     the log file.
    !
    !     Called by:     program MEAMfit
    !     Calls:         -
    !     Arguments:     ndatapoints,truedata,bestfitdata
    !     Returns:       -
    !     Files read:    -
    !     Files written: -
    !
    !     Andy Duff, Feb 2008.
    !
    !---------------------------------------------------------------------c

    use m_datapoints
    use m_optimization
    use m_geometry
    use m_filenames
    use m_atomproperties

    implicit none

    logical extraOutput
    integer idatapoint,i,j,previous,outch

    idatapoint=1
    previous=0
    do i=1,nstruct
        if (optimizeforce(i).eq.0) then
            if (extraOutput.eqv..false.) then
                if (previous.ne.1) then
                    write(outch,*)
                    write(outch,*) 'Energies:'
                    write(outch,*)
                    write(outch,*) 'Structure          fitdata', &
                        '       truedata'
                    write(outch,*) '----------------------------------------', &
                        '----'
                    previous=1
                endif
                if (useRef) then
                   write(outch,'(A25,F25.18,A2,F25.18)') strucnames(i), &
                       (fitdata(idatapoint)-fitdata(1)),'  ',(truedata(idatapoint)-truedata(1))
                else
                   write(outch,'(A25,F25.18,A2,F25.18)') strucnames(i), &
                       fitdata(idatapoint),' ',truedata(idatapoint)
                endif
            else
                if (previous.ne.1) then
                    write(outch,'(A16,A59,A57,A28)') '% id     fitdata', &
               '            truedata            pair-pot part              ', &
                        'emb-func part            smallest sepn      ..per species', &
                        ' (e.g. (1,1); (1,2); (2,2) )'
                    previous=1
                endif
                if (useRef) then
                   write(outch, & 
                '(I5,A3,F17.12,A2,F17.12,A2,F22.15,A2,F25.15,A2,F10.5,A2)',advance="no") &
                       i,'   ',(fitdata(idatapoint)-fitdata(1)),'  ',(truedata(idatapoint)-truedata(1)), &
                       '  ',sumppStr(idatapoint),'  ',summeamfStr(idatapoint),'  ', &
                       smallestsepnStr(i),'  '
                else
                   write(outch, &
                '(I5,A3,F17.12,A2,F17.12,A2,F22.15,A2,F25.15,A2,F10.5,A2)',advance="no") &
                       i,'   ',fitdata(idatapoint),' ',truedata(idatapoint), &
                       '  ',sumppStr(idatapoint),'  ',summeamfStr(idatapoint),' ', &
                       smallestsepnStr(i),'  '
                endif
                do j=1,smlsepnperspcStr_numEntries
                    write(outch,'(A1,F10.5)',advance="no") ' ',smlsepnperspcStr(j,i)
                enddo
                write(outch,*)
            endif
            idatapoint=idatapoint+1
        else
            idatapoint=idatapoint+3*gn_forces(i)
        endif
    enddo

    idatapoint=1
    previous=0
    do i=1,nstruct
        if (optimizeforce(i).eq.0) then
            idatapoint=idatapoint+1
        else
            open(78,file='forces') !save forces for use in vasp_langevin
            do j=1,gn_forces(i)
                if (optforce(j,i).eqv..true.) then
                    if (previous.eq.0) then
                        write(outch,*)
                        write(outch,*) 'Forces:'
                        write(outch,*)
                        write(outch,*) &
                            'Structure          fitdata              ', &
                            '          truedata'
                        write(outch,*) &
                            '----------------------------------------', &
                            '--------------------------------------'
                        previous=1
                    endif
                    if (previous.eq.1) then
                        !                print *,idatapoint,':'
                        write(outch,'(A25,F10.5,F10.5,F10.5,F10.5,F10.5,F10.5)') &
                            strucnames(i), &
                            fitdata(idatapoint),fitdata(idatapoint+1),fitdata(idatapoint+2), &
                            truedata(idatapoint),truedata(idatapoint+1),truedata(idatapoint+2)
                        previous=2
                    else
                        !                print *,idatapoint,':'
                        write(outch,'(A25,F10.5,F10.5,F10.5,F10.5,F10.5,F10.5)') &
                            '                    ', &
                            fitdata(idatapoint),fitdata(idatapoint+1),fitdata(idatapoint+2), &
                            truedata(idatapoint),truedata(idatapoint+1),truedata(idatapoint+2)
                    endif
                    !               fitdata(idatapoint)=-4.44089e-16
                    !               fitdata(idatapoint+1)=-4.71845e-16
                    !               fitdata(idatapoint+2)=-4.85723e-16
                    write(78,'(E14.6,E14.6,E14.6)') &
                        fitdata(idatapoint),fitdata(idatapoint+1),fitdata(idatapoint+2)
                endif
                idatapoint=idatapoint+3
            enddo
            close(78)
        endif
        if (previous.eq.2) then
            previous=1
        endif
    enddo

end subroutine writemeamdatapoints
subroutine writemeamp

    !---------------------------------------------------------------c
    !
    !     Writes MEAM parameters from the 'p' array into ouput
    !     channel 2 (should be defined w.r.t a filename prior to
    !     call)
    !
    !     Called by:     optimizeparameters
    !     Calls:         -
    !     Arguments:     maxspecies,lmax,workparameterfile,p
    !     Returns:       -
    !     Files read:    -
    !     Files written: workparameterfile
    !
    !---------------------------------------------------------------c

    use m_filenames
    use m_atomproperties
    use m_generalinfo
    use m_optimization
    use m_meamparameters

    implicit none

    integer i,j,paracounter

    !First check that non of the values are NaN
    do i=1,2*m3+(5+lm1)*m1+34*m2+12*lm1*m2
        if ((p(i).le.0d0).or.(p(i).ge.0d0)) then
            !ok...
        else
            print *,'NaN detected -not going to write', &
                'work_meam_parameters',p(i)
            !stop
        endif
    enddo

    write(2,*) lmax,' # lmax'
    write(2,*) p(1:m3),' # cmin' !cmin(maxspecies,maxspecies,maxspecies)
    write(2,*) p(1+m3:2*m3),' # cmax'    !cmax(maxspecies,maxspecies,maxspecies)
    write(2,*) envdepmeamt,' # environ dep meamts? ..followed by meamtaus:'
    write(2,*) p(2*m3+1:2*m3+lm1*m1) !meamtau(0:lmax,maxspecies)
    write(2,*) thiaccptindepndt,' # acceptor independant thi?'
    paracounter=2*m3+lm1*m1
    do i=1,maxspecies
        do j=1,maxspecies
            if ((i.eq.1).or.(thiaccptindepndt.eqv..false.)) then
                write(2,*) typethi(i,j),' #type of thi, followed by thi spec:'
                write(2,*) p(paracounter+1:paracounter+12*lm1)
            endif
            paracounter=paracounter+12*lm1
        enddo
    enddo
    write(2,*) embfunctype,' # embedding function type, followed by spec:'
    write(2,*) p(paracounter+1:paracounter+m1) !meame0(maxspecies)
    write(2,*) p(paracounter+m1+1:paracounter+2*m1)
    !meamrho0(maxspecies)
    write(2,*) p(paracounter+2*m1+1:paracounter+3*m1) !meamemb3
    write(2,*) p(paracounter+3*m1+1:paracounter+4*m1) !meamemb4
    paracounter=paracounter+4*m1
    do i=1,maxspecies
        do j=1,maxspecies
            if (j.ge.i) then
                write(2,*) typepairpot(i,j),' # pair-potential type, followed',&
                               'by spec:'
                write(2,*) p(paracounter+1:paracounter+32)
            endif
            paracounter=paracounter+32
        enddo
    enddo
    write(2,*) p(paracounter+1:paracounter+m2),' # rs' !rs(maxspecies,
    !maxspecies)
    write(2,*) p(paracounter+m2+1:paracounter+2*m2),' # rc' !rc(maxspecies,
    !maxspecies)
    write(2,*) p(paracounter+2*m2+1:paracounter+2*m2+m1),' enconst (fixed offset per atom)' !enconst
    !(maxspecies)

    if ((paracounter+2*m2+m1.ne.362).and.(m1.eq.2).and.(lmax.eq.3)) then
        print *,'Error in writemeamp subroutine regarding array limits'
        print *,'number of parameters in p=',paracounter+2*m2+m1, &
            ' not equal 362'
        stop
    endif

end subroutine writemeamp

subroutine writeoptparasSub

    !-------------------------------------------------------------------c
    !
    !     Read in the which potential parameters are to be optimized.
    !
    !     Called by:     readsettings
    !     Calls:         -
    !     Arguments:     -
    !     Returns:       -
    !     Files read:    -
    !     Files written: -
    !
    !     Andrew Duff 2014
    !
    !-------------------------------------------------------------------c

    use m_optimization
    use m_atomproperties
    use m_meamparameters

    implicit none

    integer i,isp,jsp,spcCnt,l

    do i=1,np
       if (p(i).ne.0d0) then
          freep(i)=1
       endif
    enddo

    print *
    print *,'-------------------------------------'
    print *

    write(*,'(A12)') 'PARASTOOPT={'
    do i=1,m3
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    do i=1,m3
       write(*,'(I1,A1)',advance='no') freep(i+m3),' '
    enddo
    write(*,*)
    do i=2*m3+1,2*m3+lm1*m1
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    spcCnt=0
    do isp=1,m1
       do jsp=1,m1
          do l=0,lmax
             do i=1,12
                write(*,'(I1,A1)',advance='no') &
          freep(2*m3+lm1*m1+12*(spcCnt*lm1+l)+i),' '
             enddo
             write(*,*)
          enddo
          spcCnt=spcCnt+1
       enddo
    enddo
    do i=2*m3+lm1*m1+12*lm1*m2+1,2*m3+m1 &
        +lm1*m1+12*lm1*m2
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    do i=2*m3+m1+lm1*m1+12*lm1*m2+1, &
        2*m3+2*m1+lm1*m1+12*lm1*m2
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    do i=2*m3+2*m1+lm1*m1+12*lm1*m2+1, &
        2*m3+3*m1+lm1*m1+12*lm1*m2
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    do i=2*m3+3*m1+lm1*m1+12*lm1*m2+1, &
        2*m3+4*m1+lm1*m1+12*lm1*m2
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    spccnt=0
    do isp=1,m1
       do jsp=1,m1
          do i=1,32
             write(*,'(I1,A1)',advance='no') freep(2*m3+(4+lm1)*m1+ &
                   12*lm1*m2+32*spccnt+i),' '
          enddo
          write(*,*)
          spccnt=spccnt+1
       enddo
    enddo
    !pairpotparameter(16,maxspecies,maxspecies)
    do i=m4+1,m4+m2
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    do i=m4+m2+1,m4+2*m2
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    do i=m4+2*m2+1,m4+2*m2+m1
       write(*,'(I1,A1)',advance='no') freep(i),' '
    enddo
    write(*,*)
    write(*,'(A1)') '}'
    write(*,*)
    close(1)


 !  write(*,*) freep(1:m3)  !cmin(maxspecies,maxspecies,maxspecies)
 !  write(*,*) freep(1+m3:2*m3) !cmax(maxspecies,maxspecies,maxspecies)
 !  write(*,*) freep(2*m3+1:2*m3+lm1*m1) !meamtau(0:lmax,maxspecies)
 !  write(*,*) freep(2*m3+lm1*m1+1:2*m3+lm1*m1+12*lm1*m2)
 !  !meamrhodecay(6,0:lmax,maxspecies,maxspecies)
 !  !meamrhodecay(6,0:lmax,maxspecies,maxspecies)
 !  write(*,*) freep(2*m3+lm1*m1+12*lm1*m2+1:2*m3+m1 &
 !      +lm1*m1+12*lm1*m2)   !meame0(maxspecies)
 !  write(*,*) freep(2*m3+m1+lm1*m1+12*lm1*m2+1:2*m3+ &
 !      2*m1+lm1*m1+12*lm1*m2) !meamrho0(maxspecies)
 !  write(*,*) freep(2*m3+2*m1+lm1*m1+12*lm1*m2+1:2*m3+ &
 !      3*m1+lm1*m1+12*lm1*m2) !meamemb3
 !  write(*,*) freep(2*m3+3*m1+lm1*m1+12*lm1*m2+1:2*m3+ &
 !      4*m1+lm1*m1+12*lm1*m2) !meamemb4
 !  write(*,*) freep(2*m3+(4+lm1)*m1+12*lm1*m2+1:m4)
 !  !pairpotparameter(16,maxspecies,maxspecies)
 !  write(*,*) freep(m4+1:m4+m2) !rs(maxspecies,maxspecies)
 !  write(*,*) freep(m4+m2+1:m4+2*m2) !rc(maxspecies,maxspecies)
 !  write(*,*) freep(m4+2*m2+1:m4+2*m2+m1) !enconst(maxspecies)

    print *
    print *,'Please copy and paste the above into the settings file.'
    print *,'0=no opt; 1=opt but do not randomly seed;'
    print *,'2=opt and randomly seed; 3=no opt but randomly seed'
    print *,'Ordering corresponds to that in the potparas_best files'

end subroutine writeoptparasSub
