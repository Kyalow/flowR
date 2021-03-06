% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection_module.R
\name{selectionInput}
\alias{selectionInput}
\alias{selection}
\title{Select samples and subsets from a GatingSet}
\usage{
selectionInput(id)

selection(
  input,
  output,
  session,
  rval,
  params = reactiveValues(),
  multiple_subset = TRUE
)
}
\arguments{
\item{id}{shiny id}

\item{input}{shiny input}

\item{output}{shiny output}

\item{session}{shiny session}

\item{rval}{reactivevalues object with the following elements :
\describe{
  \item{gating_set}{: a GatingSet object}
  \describe{
     \item{gate}{: a filter object (as build using polygonGate() for instance)}
     \item{parent}{: the name of the parent gate}
     }
}}

\item{params}{reactivevalues object used to initialize selected samples 
and subsets with elements (not mandatory) :
\describe{
  \item{sample}{: initially selected samples}
  \item{subset}{: initially selected subsets}
 }}

\item{multiple_subset}{logical; allow selection of multiple subsets}
}
\value{
a reactivevalues object with input values amongst which:
\describe{
  \item{sample}{: selected samples}
  \item{subset}{: selected subsets}
}
}
\description{
Select samples and subsets from a GatingSet

selection module server function
}
\examples{
\dontrun{
library(shiny)
if (interactive()){
  
  ui <- fluidPage(
    selectionInput("selection_module")
  )
  
  server <- function(input, output, session) {
    
    rval <- reactiveValues()
    params <- reactiveValues()
    
    observe({
      gs <- load_gs("./inst/ext/gs")
      rval$gating_set <- gs
      params$sample <- pData(gs)$name[2]
    })
    
    callModule(selection, "selection_module", rval, params = params)
  }
  
  shinyApp(ui, server)
  
}
}
}
