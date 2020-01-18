\name{phyloch-news}
\alias{phyloch-news} 
\title{News and changes to the phyloclim package}

\section{VERSION 0.9.8}{
  \subsection{IMPROVED FEATURES}{
    \itemize{
      \item Random presence points were drawn by \code{bg.similarity.test()} over the entire space
      defined by argument \code{env}, i.e. the bioclimatic predictors, which in most cases is 
      maybe not a reasonable choice. Now you can use the new argument \code{study.area.y = "mcp"} 
      to sample random presence points of the second species (species Y of Warren et al. 2008) 
      from minimum spanning polygon laid around the observed presence points of Y. Thanks to Frank
      Burbrink for bringing this issue up!
      }
  }
}
\section{VERSION 0.9.7}{
  \itemize{
    \item Only minor changes.
  }
}
\section{VERSION 0.9.6}{
  \itemize{
    \item Only minor changes.
  }
}
\section{VERSION 0.9.5}{
  \itemize{
    \item Only minor changes.
  }
}
\section{VERSION 0.9.4}{
  \itemize{
    \item Only minor changes.
  }
}
\section{VERSION 0.9-3}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Due to an unconsistent behavior of sapply some data 
      produced meaningless p-values in \code{age.range.correlation()}. 
      It has been replaced by lapply to be robust. 
      (thanks to Dan Warren for bug report)
    }
  }
}
\section{VERSION 0.9-2}{
  \subsection{BUG FIXES}{
    \itemize{
      \item There was a mistake in the calculation of p-values in
      \code{niche.equivalency.test()} and the print method for both
      \code{niche.quivalency.test()} and \code{bg.similarity.test()} stated the null
      hypothesis instead of the alternative hyposthesis. (thanks to 
      Guarino Colli for bug report)
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item The confidence level for CIs in \code{bg.similarity.test()} can now be
      set via the option \code{conf.level}.
      \item \code{niche.equivalency.test()} and \code{bg.similarity.test()} are now able to use
      categorical environmental covariates.
    }
  }
}
\section{VERSION 0.9-1}{
  \subsection{IMPROVED FEATURES}{
    \itemize{
      \item Calculation of nonparametric confidence intervals in 
      \code{bg.similarity.test()}.
      }
  }
}
\section{VERSION 0.9-0}{
  \subsection{IMPROVED FEATURES}{
    \itemize{
      \item The code of functions \code{niche.equivalency.test()} and \code{bg.similarity.test()}
      has been completely rewritten to match exactly the definition in 
      Warren et al. (2011). Previous versions had relied entirely on the 
      SWD-format (see MAXENT tutorial), which was quicker and easy to im-
      plement at the cost of altering the definition of the back-
      ground similarity test.
      \item \code{niche.equivalency.test()} and \code{bg.similarity.test()} now produce output
      objects of class \code{'ntest'}. For class \code{'ntest'} a print and a plot
      method have been included, the latter replacing the former 
      non-generic functions \code{plotNicheEquivalency()} and 
      \code{plotBackgroundSimilarity()}.
    }
  }
}
\section{VERSION 0.8-5}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Plotting of tiplabels in function \code{plotAncClim()} has been completey
      reworked and is now robust regarding the scale of the phylogenetic 
      tree. Have a look at \code{?plotAncClim} and note the new arguments \code{tipmode}, 
      \code{tipspace}, \code{cladespace}, and \code{nchar}. The functionality of arguments 
      respectively. Note that cladespace is still in a somewhat 
      \code{show.tips} and \code{xspace} has been included in \code{tipmode} and 
      \code{cladespace}, experimental state.
    }
  }
}
\section{VERSION 0.8-4}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Code of function \code{pno()} was reworked and is now ca. 20 percent faster
      than in version 0.8.3. This is due to the replacement of the 
      dependency on \strong{adehabitat} for raster map import by pno's own import 
      device for raster layers.
      \item A check for correct path_model argument is included in function \code{pno()}.
    }
  }
}
\section{VERSION 0.8-3}{
  \subsection{IMPROVED FEATURES}{
    \itemize{
      \item Improved documentation of \code{plotAncClim()}.
    }
  }
}
\section{VERSION 0.8-2}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Some users obtained PNOs containing NA values because the 
      distibution of NODATA values (e.g., -9999) in their bioclimatic data 
      and in their ENM models was different. This problem might arise when
      raster maps differ slightly in their coastlines. The function \code{pno()}
      will now issue a warning in these cases and treat the corresponding 
      cells as having zero probability in the ENM (e.g., the MAXENT    
      distribution).
      \item \code{plotPNO()} is now capable of handling PNO matrixes with only one 
      species.
    }
  }
}
\section{VERSION 0.8-1}{
  \itemize{
    \item Only minor changes in the documentation of PHYLOCLIM.
  }
}
\section{VERSION 0.8}{
  \subsection{BUG FIXES}{
    \itemize{
    \item The formula of Fitzpatrick and Turelli (2006) was incorrectly 
    implemented in the function \code{age.range.correlation()}. First, it 
    used 'n' (instead of n - 1) to weight the overlap statistic 'o'. 
    Second, the function was not summing 'o' correctly over all 
    species in both clades. (thanks to Pedro Abellán and 
    Ben Fitzpatrick for bug report and suggestions for a fix)
    }
  }
}
\section{VERSION 0.6}{
  \subsection{BUG FIXES}{
    \itemize{
      \item There was mistake in the formala for the calculation of 
      Hellingers distance in the original publication of Warren,
      Glor and Turelli (2008, Evolution 62:2868-2883; see also
      http://enmtools.blogspot.com/2010/09/error-in-i-metric-but-
      dont-panic.html). The corrected formula is now implemented
      in \code{getDI()}. As stressed by Dan Warren the error is only a
      scaling issue and does not affect the p-values of the tests 
      based upon it. (thanks to Pascal Title for bug report)
    }
  }
}
\section{VERSION 0.5}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Function \code{bg.similarity.test()} was mislead by different 
      choices of column names in the SWD files, e.g. "species" 
      instead of "spec". The function can now handle any column 
      name. (thanks to Pascal Title for bug report)
    }
  }
}
\section{VERSION 0.4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The functions \code{niche.equivalency.test()} and \code{bg.similarity.test()} have
      a new argument \code{mx}, which controls the number of memory that is given
      to JAVA during the calls to MAXENT.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item The functions \code{niche.equivalency.test()} and \code{bg.similarity.test()} were
      not fully Window compatible. (thanks to Alastirs Potts for the fix)
    }
  }
}
\section{VERSION 0.3}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The function \code{plotPNO()} has now an argument \code{legend.pos} to control
      the position of the legend. Alternatively, the plotting of
      the legend can now be suppressed.
    }
  }
}