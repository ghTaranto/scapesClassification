<!-- README.md is generated from README.Rmd. Please edit that file -->

# What is _scapesClassification_ ?

[_scapesClassification_](https://github.com/ghTaranto/scapesClassification "GitHub repository")
is an R-package for the classification of seascapes, landscapes and
other geo-spaces based on user-defined conditions. It allows users to
translate task-oriented views of spaces and of geographic objects into
computer representations.

# How does _scapesClassification_ work?

Geo-spaces are classified using sets of user-defined conditions. These
sets of conditions, presented in the form of [conditional
statements](<https://en.wikipedia.org/wiki/Conditional_(computer_programming)> "Wikipedia definition"),
can be applied simultaneously or sequentially. The second alternative is
generally preferred. In fact, the rationale that guided the design of
_scapesClassification_ is better implemented considering multi-step
classification processes:

1.  Geographic objects and segments of space can be identified by a
    unique class (or set of classes);

2.  Generally, it is easier to identify attributes that are distinctive
    of a portion of a class rather then identify its full range of
    attributes at once;

3.  Thus, distinctive attributes can be used to map an initial set of
    locations to a class. These locations are hereafter referred to as
    **_anchor locations_**;

4.  Then, it becomes possible to map a location to a class not only
    considering its intrinsic attributes, but also its spatial
    relationships with anchor locations;

5.  In particular, _class contiguity_ and _class continuity_ can be
    considered;

6.  **_Class contiguity._** Qualitative and/or quantitative conditions
    define the membership of a location to a class only if they occur in
    positions considered adjacent to a specific class;

7.  **_Class continuity._** When two adjacent locations are assigned to
    the same class using contiguity conditions, these conditions can be
    re-applied at positions neighboring the newly classified locations;

8.  Distinct classes can be identified repeating the above process and
    considering the relationships expected to exist among classified and
    unclassified locations.
