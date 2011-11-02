\name{SDiscReportHead}
\alias{SDiscReportHead}
\title{SDiscReportHead}
\description{LaTeX header for Sweave reporting.}
\usage{SDiscReportHead(LO='Report', packages=list('Sweave','amsmath', 'underscore', 'setspace', 'pdflscape',
   'multirow', 'glossaries'=c('toc', 'acronym', 'xindy'), 'multicol', inputenc='latin1', babel='english',
   'pdfpages', caption=c('small','bf'), 'graphicx','fancyhdr', 'lastpage','longtable','color',xcolor='table',
   geometry=c('left=1.25cm','top=2cm','right=1.25cm','bottom=2cm'), hyperref=c('colorlinks=true',
   'citecolor=blueDoc', 'filecolor=blueDoc', 'linkcolor=blueDoc',
   'urlcolor=blueDoc'),'lscape','sectsty','colortbl','wrapfig','array'),
   author=list(name='MyName',email='MyEmail',address='address'))}
\arguments{
   \item{LO}{ the top left page header (default is "Report")}
   \item{packages}{ named list of packages. List names refer to LaTeX packages, and are passed to the
   \code{usepackage} command. For each name, the character vector defines the list of options (squared brackets) of
   the LaTeX package, comma separated. }
   \item{author}{ defines the author list; for one author, a list of three elements (\code{name}, \code{email},
      \code{address}); for more authors, a list element of lists composed of (\code{name}, \code{email},
      \code{address}).}
}
\author{Fabrice Colas}
\keyword{print}
\alias{SDiscReportHead}
