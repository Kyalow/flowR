% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flowR.R
\name{add_columns_from_metadata}
\alias{add_columns_from_metadata}
\title{Add metadata columns}
\usage{
add_columns_from_metadata(df, metadata)
}
\arguments{
\item{df}{data.frame with a column \code{name} used to map metadata.
Metadata should also contain a column \code{name}}

\item{metadata}{a data.frame containing metadata associated to samples. 
Must have a column \code{name} used for mapping.}
}
\value{
a data.frame with additional columns
}
\description{
Add metadata columns
}
