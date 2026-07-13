# === HD 2×2 PANELS FROM EXISTING PNGs (MISSPEC) ===============================
# Uses the PNGs already in complete_misspec_results to build four big panels:
#   1) full • exposed_to_infected
#   2) full • susceptible_to_exposed
#   3) partial • exposed_to_infected
#   4) partial • susceptible_to_exposed

# install.packages("magick") # if needed
library(magick)

# ── CONFIG: point to your misspec folder ────────────────────────────────────────
base_dir       <- "~/Rt_calculation_serial_interval/complete_misspec_results/"
out_dir        <- file.path(base_dir, "grids_hd")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# visual sizing (tweak if desired)
per_plot_width <- 4300   # width of each source PNG in the grid (px)
frame_px       <- 50     # white frame around each cell (px)
gap_px         <- 50     # vertical space between top/bottom rows (px)
title_height   <- 160    # title bar height (px)
title_size     <- 64     # title font size
label_w        <- 820    # R0 label box width (px)
label_h        <- 160    # R0 label box height (px)
label_size     <- 66     # R0 label font size
label_bg       <- "rgba(255,255,255,0.90)"
bg_color       <- "white"

# order in the 2×2 grid: (1.5 | 2) / (3 | 5)
r0_vals <- c("1.5","2","3","5")

# ── HELPERS ─────────────────────────────────────────────────────────────────────
find_png <- function(model_type, r0, transition) {
  candidates <- c(
    file.path(base_dir, sprintf("combined_plot_%s_R0_%s_%s.png", model_type, r0, transition)),
    file.path(base_dir, sprintf("combined_plot_%s_R0_%s_%s.png", model_type, gsub("\\.", "_", r0), transition))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

prep_cell <- function(path, r0_text) {
  im <- image_read(path)
  im <- image_scale(im, paste0(per_plot_width, "x"))
  im <- image_border(im, color = bg_color, geometry = sprintf("%dx%d", frame_px, frame_px))
  
  # R0 label overlay (top-left)
  label <- image_blank(width = label_w, height = label_h, color = label_bg)
  label <- image_annotate(label, paste0("R0 = ", r0_text),
                          gravity = "center", size = label_size, color = "black")
  image_composite(im, label, offset = sprintf("+%d+%d", frame_px + 20, frame_px + 20))
}

make_grid <- function(model_type, transition, title) {
  files <- vapply(r0_vals, function(rv) find_png(model_type, rv, transition), character(1))
  if (anyNA(files)) {
    stop("Missing PNG(s) for ", model_type, " / ", transition,
         " at R0: ", paste(r0_vals[is.na(files)], collapse = ", "),
         "\nLooked in: ", normalizePath(base_dir))
  }
  
  cells <- Map(prep_cell, files, r0_vals)
  
  row1 <- image_append(image_join(cells[[1]], cells[[2]]))   # left-right
  row2 <- image_append(image_join(cells[[3]], cells[[4]]))
  
  w <- max(image_info(row1)$width, image_info(row2)$width)
  spacer <- image_blank(width = w, height = gap_px, color = bg_color)
  grid <- image_append(image_join(row1, spacer, row2), stack = TRUE)
  
  title_bar <- image_blank(width = image_info(grid)$width, height = title_height, color = bg_color)
  title_bar <- image_annotate(title_bar, title, gravity = "center", size = title_size, color = "black")
  grid <- image_append(image_join(title_bar, grid), stack = TRUE)
  
  out_file <- file.path(out_dir, sprintf("grid_%s_%s_hd.png", model_type, transition))
  image_write(grid, out_file, format = "png")
  cat("✅ Saved:", out_file, " (", image_info(grid)$width, "×", image_info(grid)$height, "px)\n", sep = "")
  out_file
}

# ── BUILD THE FOUR MISSPEC PANELS ───────────────────────────────────────────────
g1 <- make_grid("full",    "exposed_to_infected",    "FULL • Exposed → Infected")
g2 <- make_grid("full",    "susceptible_to_exposed", "FULL • Susceptible → Exposed")
g3 <- make_grid("partial", "exposed_to_infected",    "PARTIAL • Exposed → Infected")
g4 <- make_grid("partial", "susceptible_to_exposed", "PARTIAL • Susceptible → Exposed")

c(g1, g2, g3, g4)
