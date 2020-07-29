DCM\_investigation
================
MAST
February 20, 2020

-   [Evolution of the Black Sea Deep Chlorophyll Maximum](#evolution-of-the-black-sea-deep-chlorophyll-maximum)
-   [Seasonal Time-Series](#seasonal-time-series)
-   [Back Scattering Validation of the DCM](#back-scattering-validation-of-the-dcm)
-   [Depth horizons](#depth-horizons)
    -   [Density Horizons](#density-horizons)
    -   [Vertical distribution of the CHL content](#vertical-distribution-of-the-chl-content)
    -   [Concentrations](#concentrations)
    -   [Navarro's Ratio](#navarros-ratio)
    -   [and what about light at DCM ?](#and-what-about-light-at-dcm)
-   [Potential Additional sections.](#potential-additional-sections.)
    -   [Spatial](#spatial)
    -   [Inter-annual variability.](#inter-annual-variability.)

Evolution of the Black Sea Deep Chlorophyll Maximum
===================================================

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig3-1.png)

Seasonal Time-Series
====================

Back Scattering Validation of the DCM
=====================================

    ## Loading required package: reshape2

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig9-1.png)

Depth horizons
==============

    ## notch went outside hinges. Try setting notch=FALSE.

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig5-1.png)

Density Horizons
----------------

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig6-1.png)

Vertical distribution of the CHL content
----------------------------------------

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig8-1.png)

Concentrations
--------------

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig7-1.png)

Navarro's Ratio
---------------

``` r
require(scales)

MLDdf_info2 <- ddply(tmp_fitco, .(year,platform), summarize,
                maxMLD = max(MLD, na.rm = T),
                sigmaMaxMLD = max(sigmaMLD, na.rm = T))

sigma_ratio2 <- ddply(tmp_fitco, .(platform, juld), mutate,
   sigmaMaxMLD = MLDdf_info2$sigmaMaxMLD[which((MLDdf_info2$year == unique(year))
                                              & (MLDdf_info2$platform == unique(platform)))],
   ratio = sigmaMAX/sigmaMaxMLD)
```

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/ratio_spatial2b-1.png)

and what about light at DCM ?
-----------------------------

    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.
    ## notch went outside hinges. Try setting notch=FALSE.

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/PARatDCM2BIS-1.png)

Potential Additional sections.
==============================

Spatial
-------

    ## Loading required package: ggmap

    ## Google's Terms of Service: https://cloud.google.com/maps-platform/terms/.

    ## Please cite ggmap if you use it! See citation("ggmap") for details.

    ## Source : http://tile.stamen.com/terrain/6/36/22.png

    ## Source : http://tile.stamen.com/terrain/6/37/22.png

    ## Source : http://tile.stamen.com/terrain/6/38/22.png

    ## Source : http://tile.stamen.com/terrain/6/39/22.png

    ## Source : http://tile.stamen.com/terrain/6/36/23.png

    ## Source : http://tile.stamen.com/terrain/6/37/23.png

    ## Source : http://tile.stamen.com/terrain/6/38/23.png

    ## Source : http://tile.stamen.com/terrain/6/39/23.png

    ## Source : http://tile.stamen.com/terrain/6/36/24.png

    ## Source : http://tile.stamen.com/terrain/6/37/24.png

    ## Source : http://tile.stamen.com/terrain/6/38/24.png

    ## Source : http://tile.stamen.com/terrain/6/39/24.png

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig4-1.png)

Inter-annual variability.
-------------------------

Here we repeat all plots made above, but with single lines marking interannual variability.

![](article_figure_ARTHUR_CAPET_files/figure-markdown_github/Fig11-1.png)
