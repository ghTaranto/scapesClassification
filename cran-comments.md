Dear CRAN maintainers,

I am submitting a new package called 'scapesClassification' 

Thanks you,

Gerald H. Taranto

## Release summary

* This is a new release.

## Test environments

_devtools::check()_

* Local x86_64-w64-mingw32, R 4.0.5

_devtools::check_win_release(), devtools::check_win_devel() & 
devtools::check_mac_release()_

* x86_64-w64-mingw32 64-bit, R 4.1.3

* x86_64-w64-mingw32 (64-bit), R-devel

* aarch64-apple-darwin20 (64-bit), R 4.1.1

_rhub::check_for_cran()_

* Windows Server 2022, R-devel, 64 bit

* Ubuntu Linux 20.04.1 LTS, R-release, GCC

* Fedora Linux, R-devel, clang, gfortran

_rhub::check_on_debian()_

* Debian Linux, R-release, GCC

## R CMD check results
  
There was 1 warning common to many environments

> checking for executable files ... WARNING
  Found the following executable file:
    inst/extdata/Azores.dbf

According to dbf specification (found in this link https://www.dbase.com/Knowledgebase/INT/db7_file_fmt.htm ) a DBASE level 5 file, last updated in 2022 (makes second byte 122, that's 122 years after 1900 :P ), matches the above 2-byte signature, so it gets misidentified as executable (link)[https://mac.r-project.org/macbuilder/results/1647016745-f231131c578998f4/]
