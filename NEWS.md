# dispersionIndicators 0.2.0
* New functionnalities for plotting convex hulls
  * Display of convex hulls indicators in scientifc notations on plots
  * Batch names are displayed on shared plots
  * To avoid legend overlapping on the plot, on shared plots, the legend is displayed on a separated pdf page
  * Possibility to choose a compact order on the x-axis to better visualise batches. Available for shared plot and plot 
    all in global mode.
* Bug of integer overflow fixed


# dispersionIndicators 0.1.6
* New convex plot function: `shared_single_plot_convex_hulls`
  * Allows for the creation of a single plot displaying convex hulls for multiple variables.

# dispersionIndicators 0.1.5
* Added URL to description

# dispersionIndicators 0.1.4
* Fixed some issues encountered on CRAN checking pipelines.
  * The test on output files for convex functions was disabled for unstable version of R under windows OS, due to 
    too much instabilities observed on CRAN pipelines.

# dispersionIndicators 0.1.3
* Fixed some issues encountered on CRAN checking pipelines.
  * Replace test of outfile files for convex functions from pdf files to png
    files (obtained from the pdf files). This avoids issues of unstabilities of
    pdf content (due to e.g. software versions).

# dispersionIndicators 0.1.2
* Fixed some issues pointed out by CRAN team members.
  * Make tests pass for mac OS (issue with generated files checksum)

# dispersionIndicators 0.1.1
* Fixed some issues pointed out by CRAN team members.
  * Writing functions in examples/vignettes/tests to tmp files
  * More explicit copyrights
  * Adding references in the DESCRIPTION file


# dispersionIndicators 0.1.0
* Initial package version.











