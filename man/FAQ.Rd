\name{FAQ}
\alias{FAQ}

\title{ Frequently Asked Questions, And Others }

\description{
  
  A place for hints and miscellaneous advice. 
  
}

\section{How do I install and start secr?}{

Follow the usual procedure for installing from CRAN archive (see menu
item Packages | Install package(s)... in Windows). You also need to get
the package \pkg{abind} from CRAN.

Other required packages (\pkg{MASS}, \pkg{nlme}, \pkg{stats}) should be
available as part of your \R installation.

Like other contributed packages, \bold{secr} needs to be loaded before
each use e.g.,\code{library(secr)}.

You can learn about changes in the current version with
\code{news(package = "secr")}.

}

\section{ How can I get help?}{

There are three general ways of displaying documentation from within \R.
Firstly, you can bring up help pages for particular functions from the
command prompt.  For example:

\code{?secr}
or
\code{?secr.fit}

Secondly, help.search() lets you ask for a list of the help pages on a
vague topic (or just use ?? at the prompt).  For example:

\code{?? "linear models"}

Thirdly, you can display various \pkg{secr} documents listed in
\code{\link{secr-package}}.

Tip: to search all secr help pages open the pdf version of the manual in
Acrobat Reader (\href{https://www.otago.ac.nz/density/pdfs/secr-manual.pdf}{secr-manual.pdf}; see also ?secr) and use
<ctrl> F.

There is a support forum at \url{http://www.phidot.org/forum/} under
`DENSITY|secr' and another at \href{https://groups.google.com/forum/#!forum/secrgroup}{secrgroup}. See below for
more \R tips. Some specific problems with \code{secr.fit} are covered in
\link{Troubleshooting}.

}

\section{ How should I report a problem? }{

If you get really stuck or find something you think is a bug then please
report the problem to one of the online lists.

You may be asked to send an actual dataset - ideally, the simplest one
that exhibits the problem. Use \code{\link{save}} to wrap
several \R objects together in one .RData file, e.g.,
\code{save("captdata", "secrdemo.0", "secrdemo.b", file =
"mydata.RData")}. Also, paste into the text of your message the output
from \code{packageDescription( "secr" )}.

}

\section{ Why do I get different answers from secr and Density? }{

Strictly speaking, this should not happen if you have specified the same
model and likelihood, although you may see a little variation due to the
different maximization algorithms. Likelihoods (and estimates) may
differ if you use different integration meshes (habitat masks), which
can easily happen because the programs differ in how they set up the
mesh.  If you want to make a precise comparison, save the Density mesh
to a file and read it into \pkg{secr}, or vice versa.

Extreme data, especially rare long-distance movements, may be handled
differently by the two programs. The `minprob' component of the
`details' argument of \code{secr.fit} sets a lower threshold of
probability for capture histories (smaller values are all set to
minprob), whereas Density has no explicit limit.

}

\section{ How can I speed up model fitting and model selection?}{

  There are many ways - see \link{Speed tips} and 
  \href{https://www.otago.ac.nz/density/pdfs/secr-troubleshooting.pdf}{secr-troubleshooting.pdf}.

}

\section{ Does secr use multiple cores?}{
  
  Some computations can be run in parallel on multiple processors (most
  desktops these days have multiple cores). Likelihood calculations in 
  \code{secr.fit} assign capture histories to multiple parallel threads 
  whenever possible. 
  
}

\section{ Can a model use detector-level covariates that vary over
  time?}{

  Yes. See ?timevaryingcov. However, a more direct way to control for
  varying effort is provided - see the `\link{usage}' atribute, which
  now allows a continuous measure of effort
  (\href{https://www.otago.ac.nz/density/pdfs/secr-varyingeffort.pdf}{secr-varyingeffort.pdf}). A tip: covariate models fit
  more quickly when the covariate takes only a few different values.
  
  }

\section{Things You Might Need To Know About \R}{

The function \code{findFn} in package \pkg{sos} lets you search CRAN for
R functions by matching text in their documentation.

There is now a vast amount of \R advice available on the web. For the
terminally frustrated, `R inferno' by Patrick Burns is recommended
(\url{https://www.burns-stat.com/pages/Tutor/R_inferno.pdf}). "If you are using
R and you think you're in hell, this is a map for you".

Method functions for S3 classes cannot be listed in the usual way by
typing the function name at the \R prompt because they are `hidden' in a
namespace. Get around this with getAnywhere(). For example:

\code{getAnywhere(print.secr)}

R objects have `attributes' that usually are kept out of sight.
Important attributes are `class' (all objects), `dim' (matrices and
arrays) and `names' (lists). \pkg{secr} hides quite a lot of useful data
as named `attributes'. Usually you will use summary and extraction
methods (\code{traps}, \code{covariates}, \code{usage} etc.) to view and change
the attributes of the various classes of object in \pkg{secr}. If you're
curious, you can reveal the lot with `attributes'.  For example, with
the demonstration capture history data `captdata':


\code{traps(captdata)       ## extraction method for `traps'}

\code{attributes(captdata)  ## all attributes}

Also, the function \code{str} provides a compact summary of any object:

\code{str(captdata)}

}

\references{
  
Claeskens, G. and Hjort N. L. (2008) \emph{Model Selection and Model
Averaging}. Cambridge: Cambridge University Press.

}

\keyword{ manip }

