% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flowR.R
\name{add_gates_flowCore}
\alias{add_gates_flowCore}
\title{Add gates from a gating hierarchy to a GatingSet}
\usage{
add_gates_flowCore(gs, gates)
}
\arguments{
\item{gs}{a GatingSet}

\item{gates}{a named list representing the gating hierarchy. 
Each element must have a field 'gate' with a flowCore filter object 
and a field 'parent' with the name of its parent gate.}
}
\description{
Add gates from a gating hierarchy to a GatingSet
}
