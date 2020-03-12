library(flowCore)
library(premessa)

# read barcode (needed)
# data("GvHD")
# 
# beads <- premessa::read_barcode_key("/mnt/NAS7/Workspace/hammamiy/data_premasse/Fluidigm_20plex_barcode_key.csv")
# file <- flowCore::read.FCS("/mnt/NAS7/Workspace/hammamiy/data_premasse/20120222_cells_found.fcs")
# 
# debarcoding <- premessa::debarcode_data(m = a@exprs, bc.key = beads)
# 
# premessa::plot_barcode_separation(debarcoding)
# premessa::plot_separation_histogram(debarcoding)
# premessa::plot_barcode_channels_intensities(exprs(a), debarcoding$bc.channels)
# premessa::plot_barcode_yields(debarcoding, 0.01)
# 
# sep <- premessa::get_assignments_at_threshold(debarcoding, sep.threshold = 0.1)
# 
# 
# mahalanobis <- premessa::get_mahalanobis_distance(a@exprs, debarcoding, 0.1)
# 
# a[, debarcoding$bc.channels]
# a@exprs[, debarcoding$bc.channels]
# plot_beads_over_time(beads.data = a@exprs, beads.normed = debarcoding$m.normed, beads.cols = debarcoding$bc.channels)
# 
# 
# 
# test <- premessa::get_assignments_at_threshold(debarcoding, 0.1, 0.95 ,mahalanobis)
# test
# 
# premessa:::debarcoder_GUI()


################################################ 
#### debarcoding example #######################

# barcode + fcs 
bead_sample <- premessa::read_barcode_key("/mnt/NAS7/Workspace/hammamiy/data_premasse/sample_barcode_key.csv")
sample <- flowCore::read.FCS("/mnt/NAS7/Workspace/hammamiy/data_premasse/sample_barcoded_file.fcs")

#save debarcoding in new file ccalled debarcoded
output.dir <- file.path("/mnt/NAS7/Workspace/hammamiy/data_premasse", "debarcoded")
dir.create(output.dir, recursive = T)
output.basename <- tools::file_path_sans_ext(basename(sample@description$FILENAME))


# select parameters 
sep.threshold <- 0.3 
label  <- row.names(bead_sample)


# result bc.res
sample_debarcoded <- premessa:::debarcode_data(sample@exprs, bead_sample) 

d_mahalanobis <- premessa:::get_mahalanobis_distance(sample@exprs, sample_debarcoded, sep.threshold = sep.threshold)
d_mahalanobis
# if separtion plot barcode separation (freq population & seuil de séparation)
premessa::plot_barcode_separation(sample_debarcoded, sep.threshold = sep.threshold)

# else  cas event
premessa::plot_barcode_yields(sample_debarcoded, sep.threshold = sep.threshold, mahal.threshold = 0.3, mahal.dist = 0.1)

mahal.dist <- cbind(sample@exprs, mahal.dist = d_mahalanobis)

sel.rows <- premessa:::get_sample_idx(label = label[1], bc.results = sample_debarcoded, sep.threshold = sep.threshold, mahal.threshold = 1, mahal.dist = 1)
m <- mahal.dist[sel.rows, ]


# (if dans le else) plot rescale event
premessa::plot_barcode_channels_intensities(sample@exprs, sample_debarcoded$bc.channels, sample_debarcoded$m.normed[sel.rows,])

# else if (single biaxial) 
x <- sample_debarcoded$bc.channels
y <- sample_debarcoded$bc.channels

premessa:::plot_color_coded_biaxial(m, x[1],
                                    y[2], "mahal.dist")

# all barcodes 

premessa::plot_all_barcode_biaxials(m, sample_debarcoded$bc.channels)


# save debarcoding 

premessa::debarcode_fcs(sample, bead_sample, output.dir = output.dir, output.basename = output.basename, sep.threshold = sep.threshold, mahal.dist.threshold = 30)



### Normalisation #########################################################################

norm <- flowCore::read.FCS("/mnt/NAS7/Workspace/hammamiy/data_premasse/20120222_cells_found.fcs")
path.remove <- "/mnt/NAS7/Workspace/hammamiy/data_premasse/20120222_cells_found.fcs"

normed.dir <- file.path("/mnt/NAS7/Workspace/hammamiy/data_premasse", "normed")
dir.create(normed.dir, recursive = T)

out.dir <- file.path("/mnt/NAS7/Workspace/hammamiy/data_premasse", "bead_removed")
dir.create(out.dir, recursive = T)
# parameter 
beadremovalui_cutoff <- 2 # entre 0 et 20
bead_type <- c("Fluidigm Beads (140,151,153,165,175)", "Beta Beads (139,141,159,169,175)")
beadremovalui.plots.number <- 3

bead.types <- premessa:::get_beads_type_from_description(bead_type)
# bead.types

# find channel 
ch_find <- premessa:::find_beads_channels_names(norm, bead.type = bead.types[2])
ch_find
# plot 
combs <- rep(ch_find, length.out = beadremovalui.plots.number * 2)

m <- flowCore::exprs(norm)

lapply(seq(1, length(combs), 2), function(i) {
  plot.idx <- ceiling(i / 2)
  plot.output <- premessa:::plot_distance_from_beads(m, combs[i], combs[i + 1])
 
})

colnames(m)
premessa:::plot_distance_from_beads(norm@exprs, ch_find[1], ch_find[5])

# remove bead from one file
premessa::remove_beads_from_file(path.remove, beadremovalui_cutoff, out.dir)


# normalization 

# visualize bead 

dna.col <- premessa:::find_dna_channel(norm)
bead.types <- premessa:::get_beads_type_from_description(bead_type)
beads.cols <- premessa:::find_bead_channels(norm, bead.types)
beads.cols.names <- premessa:::get_parameter_name(norm, beads.cols)

beads.gates<- premessa:::get_initial_beads_gates(norm)
beads.gates


sel <- premessa:::identify_beads(norm@exprs, beads.gates, beads.cols.names, dna.col)

do_plot_outputs <- function(sel.beads = NULL) {
  
  colors <- rep("black", nrow(m))
  if(!is.null(sel.beads))
    colors[sel.beads] <- "red"
  
  #Needs to be in lapply to work
  #see https://github.com/rstudio/shiny/issues/532
  lapply(1:length(beads.cols), function(i) {
    xAxisName <- premessa:::get_parameter_name(norm, beads.cols[i])
    yAxisName <- premessa:::get_parameter_name(norm, dna.col)
    
    
      list(
        x = m[, beads.cols[i]],
        y = m[, dna.col],
        color = colors,
        xAxisName = xAxisName,
        yAxisName = yAxisName,
        file = norm,
        channelGates = beads.gates$sample[[xAxisName]]
      )

  })
}

do_plot_outputs(sel)


premessa::normalize_folder(wd = "/mnt/NAS7/Workspace/hammamiy/data_premasse", output.dir.name =  "normed", beads.gates = beads.gates, beads.type = bead.types[2], baseline = NULL)

beads.gates[1]
