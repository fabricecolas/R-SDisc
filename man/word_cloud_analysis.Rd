\name{word_cloud_analysis}
\alias{word_cloud_analysis}
\title{ Association analysis between the features and a nominal attribute
shaped as a word cloud}
\description{ This method makes an HTML word cloud that can highlight the
features that most discriminate a nominal attribute, according to the chi2
test of association. }
\usage{ word_cloud_analysis(m, convert2nominal=NA, title=NA, filename="WC",
test_threshold=0.05, magn=0.05, css_style="font-family:times;",
colors=c('#fddbc7','#d1e5f0')) }
\arguments{
   \item{m}{ A data matrix,}
   \item{convert2nominal}{ The continuous and ordinal features to convert to
   nominal values with respect to the median of the feature distribution (see
   \code{to_binary}),}
   \item{title}{ The title of the word cloud,}
   \item{filename}{ The filename of the word cloud,}
   \item{test_threshold}{ The threshold above which p-values are not reported,}
   \item{magn}{ A magnitude factor that controls the span of the character
   size,}
   \item{css_style}{ The HTML-css style of each word,}
   \item{colors}{ The colors to be used.}
}
\author{ Fabrice Colas }
