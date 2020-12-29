!--------------------  M O D U L E     I N F O  ----------------------c
!
!
module m_generalinfo
    !General stuff that is carried from the main header routines
    logical readpotfile,settingsfileexist,writeExPot,contjob,writePot
    integer log10
    real(8) p_rmax !longest interatomic distance to be considered 6
    real(8), parameter:: smllnum=1.120d-16 !Smallest number for which
    ! 1 + smllnum /ne 1
    integer,parameter:: nsigdig=floor(-log10(smllnum))
    integer tstart
end module m_generalinfo
