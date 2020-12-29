module m_filenames
    !Contains the number of poscar files to be read in. Also contains
    !the names of the following files: the poscar and outcar files;
    !the files containing the input, output and workspace MEAM
    !parameters; a file which can be used to store the atomic
    !positions of one of the poscar files in xyz format; a settings
    !file; a file containing a list of all of the poscar files; a
    !file containing a list of all of the outcar files.
    logical verbose
    integer nposcarfiles
    character*80, allocatable:: poscarfiles(:)
    character*80, allocatable:: strucnames(:)
    character*80, allocatable:: outcarfiles(:)
    logical, allocatable:: vasprun(:)
    character*80 startparameterfile,crystalstrucfile,settingsfile, &
        fitdbse,lookuptablefile
end module m_filenames
