module m_optimization
    integer n_optfunc,noptfuncstore,embfuncRandType!,n_optdiffint
    integer np,nsteps         !number of meam parameters and number of
    !changes in the meam parameters allowed
    integer startparas        !1:start meam parameters taken from
    !startmeamparameters file. 2:chosen
    !randomly. 3:from genetic algorithm.
    integer nRandGen !Number of randomly generated parameters
    logical firstpraxiscall,optimizestruc,noOpt,opt_failed,cutoffopt, &
            forcefit,printoptparas,readParasOnebyone,genAlgo,fixPotIn,useRef
    logical, allocatable:: positivep(:)
    integer, allocatable:: failedtries(:),freep(:)
    integer, allocatable:: nconfig(:)
    real(8), allocatable:: p(:),scaledstep(:),p_orig(:),p_saved(:,:)
    real(8), parameter:: poptStrt=10d0
    real(8), allocatable:: fitdata(:) !size: ndatapoints
    real(8), allocatable:: bestfitdata(:) !size: ndatapoints
    real(8) avgEn,avgFrc,varEn,varFrc
    integer nEn,nFrcComp,stopTime
    real(8), allocatable:: sumppStr(:),summeamfStr(:)
    real(8), allocatable:: weights(:) !weights for the terms in the
    !optimization function. size: nstruc
    real(8), allocatable:: bestoptfuncs(:),timeoptfunc(:)
    real(8) funcforc,funcen,maxoptfuncallowed, &
        nSMAvalue,upplimoptfunc,nFX,lowestoptfuncrand,optdiff, &
        optAcc,upplimoptfunc_GA,optfunc_err
    integer maxfuncevals
    real(8) cutoffMin,cutoffMax,CutoffPenalty,CutoffPenCoeff
    real(8), parameter:: cutoffMinLim=1.5d0
    !real(8), allocatable:: pairpotStr(:,:,:)
    !integer splnNvals

    !Parameters used to randomly generate the initial MEAM parameters
    integer, allocatable:: meamtau_minorder(:,:), &
        meamtau_maxorder(:,:), &
        meamrhodecay_negvals(:,:,:), &
        meamrhodecay_minorder(:,:,:,:), &
        meamrhodecay_maxorder(:,:,:,:), &
        pairpotparameter_negvals(:,:), &
        pairpotparameter_minorder(:,:,:), &
        pairpotparameter_maxorder(:,:,:), &
        nfreeppairpot(:,:), &
        minordervaluepairpot(:),maxordervaluepairpot(:)
    real(8), allocatable:: meamrhodecay_minradius(:,:,:), &
        meamrhodecay_maxradius(:,:,:), &
        pairpotparameter_minradius(:,:), &
        pairpotparameter_maxradius(:,:), &
        meame0_minorder(:), meame0_maxorder(:), &
        enconst_minorder(:), enconst_maxorder(:), &
        meamrho0_minorder(:),meamrho0_maxorder(:), &
        meamemb3_minorder(:),meamemb3_maxorder(:), &
        meamemb4_minorder(:),meamemb4_maxorder(:)
    !Default values for setting up random parameters:
    integer, parameter:: meamtau_minorder_default = -1
    integer, parameter:: meamtau_maxorder_default = 1
    integer, parameter:: meamrhodecay_negvals_default = 3
    integer, parameter:: meamrhodecay_minorder_default = 0
    integer, parameter:: meamrhodecay_maxorder_default = 1
    integer, parameter:: meame0_minorder_default = 0
    integer, parameter:: meame0_maxorder_default = 1
    integer, parameter:: meamrho0_minorder_default = -5
    integer, parameter:: meamrho0_maxorder_default = -2
    integer, parameter:: meamemb3_minorder_default = -9
    integer, parameter:: meamemb3_maxorder_default = -6
    integer, parameter:: meamemb4_minorder_default = -1
    integer, parameter:: meamemb4_maxorder_default = 2
    integer, parameter:: pairpotparameter_negvals_default = 3
    integer, parameter:: pairpotparameter_minorder_default = 0
    integer, parameter:: pairpotparameter_maxorder_default = 1
    integer, parameter:: enconst_minorder_default = 0
    integer, parameter:: enconst_maxorder_default = 1
    logical:: readparasfromsettings

    !Genetic algo
    real(8) probRand,probPot1

end module m_optimization

