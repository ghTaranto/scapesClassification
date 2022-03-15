Dear CRAN maintainers,

I am re-submitting the package 'scapesClassification' 

after responding to your comments

Thanks you,

Gerald H. Taranto

## Re-submission

* Please do not start the description with "This package", package name, title or similar.

> Done

* References describing the methods in your package.

> At the moment there are no external references as I developed all algorithms from scatch. 
I'm working on a peer reviewd article I'll add as a reference as soon as I pass the review
process.

* Use of `\dontrun{}` in examples.

> I removed all the examples using dontrun

* Please always make sure to reset to user's options()

> Done

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
  
There was **1 WARNING** common to many environments:

* checking for executable files ... WARNING, Found the following executable file: inst/extdata/Azores.dbf

[Misidentified as executable](https://mac.r-project.org/macbuilder/results/1647016745-f231131c578998f4/):

> According to dbf specification (found in this link https://www.dbase.com/Knowledgebase/INT/db7_file_fmt.htm ) a DBASE level 5 file, last updated in 2022 (makes second byte 122, that's 122 years after 1900 :P ), matches the above 2-byte signature, so it gets misidentified as executable. **LINKS:** [link1](https://bugs.astron.com/view.php?id=316), [link2](https://stat.ethz.ch/pipermail/r-package-devel/2022q1/007722.html) and [link3](https://stackoverflow.com/questions/70713010/convincing-r-that-the-dbf-file-associated-with-a-shp-file-is-not-an-executable). 

> inst/extdata/Azores.dbf is actually a shape file I use to provide a [working example](https://ghtaranto.github.io/scapesClassification/articles/ghp/scapesClassification_02_2_ISU.html#anchor-cells).

There were **2 NOTES**:

* New submission

> No comment on my part. 

* checking for detritus in the temp directory ... NOTE, Found the following files/directories:
    'lastMiKTeXException'
    
> This note appeared only on Windows Server 2022, R-devel, 64 bit; apparently it is a bug/crash in miktex. [LINK](https://githubhelp.com/r-hub/rhub/issues/503)
