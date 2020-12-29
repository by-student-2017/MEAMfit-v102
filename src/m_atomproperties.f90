module m_atomproperties
    integer maxspecies  !The number of different species found across
    !the poscar files.
    integer m1,m2,m3,lm1,m4
    integer z2species(112) !This converts Z to species, where the
    integer speciestoz(10)
    !species are labeled as 1,2,... according
    !to which order they were found in the
    !poscar files.
    character*2 element(112)
    data element/'H','He','Li','Be','B','C','N','O','F','Ne', &
        'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti', &
        'V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se', &
        'Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd', &
        'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce', &
        'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
        'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb', &
        'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu', &
        'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg', &
        'Bh','Hs','Mt','Ds','Rg','Cn'/
    !Atomic masses in units of u (grams mol^-1). Data taken from CIAAW
    !(http://www.ciaaw.org/atomic-weights.htm#m). Where two values were 
    !given, the first is included here, and where no mass is available
    !a zero is included. (Note: these are suitable for LAMMPS when using
    !metal units)
    real(8) mass(112)
    data mass/1.00784d0,4.002602d0,6.938d0,9.0121831d0,10.806d0,12.0096d0, &
    14.00643d0,15.99903d0,18.99840316d0, 20.1797d0,22.98976928d0,24.304d0, &
    26.9815385d0,28.084d0,30.97376199d0, 32.059d0,35.446d0,39.948d0,39.0983d0,&
    40.078d0,44.955908d0,47.867d0,50.9415d0,51.9961d0,54.938044d0,55.845d0, &
    58.933194d0,58.6934d0,63.546d0,65.38d0,69.723d0,72.630d0,74.921595d0, &
    78.971d0,79.901d0,83.798d0,85.4678d0,87.62d0,88.90584d0,91.224d0,92.90637d0,&
    95.95d0, 0d0 ,101.07d0,102.90550d0,106.42d0,107.8682d0,112.414d0,114.818d0, &
    118.710d0,121.760d0,127.60d0,126.90447d0,131.293d0,132.9054519d0,137.327d0, &
    138.90547d0,140.116d0,140.90766d0,144.242d0, 0d0 ,150.36d0,151.964d0, &
    157.25d0,158.92535d0,162.500d0,164.93033d0,167.259d0,168.93422d0,173.054d0, &
    174.9668d0,178.49d0,180.94788d0,183.84d0,186.207d0,190.23d0,192.217d0, &
    195.084d0,196.966569d0,200.592d0,204.382d0,207.2d0,208.98040d0, 0d0 , 0d0, &
    0d0 , 0d0 , 0d0 , 0d0 , 232.0377d0, 231.03588d0, 238.02891, 0d0, 0d0, 0d0, &
    0d0 , 0d0 , 0d0 , 0d0 , 0d0 , 0d0 , 0d0 , 0d0 , 0d0 , 0d0 , 0d0 , 0d0 , &
    0d0 , 0d0 , 0d0 , 0d0 , 0d0/

end module m_atomproperties
