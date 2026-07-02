library(officer)
library(here)

out_path <- here("cig_filter_study_design.pptx")

prs <- read_pptx()

# ── Slide dimensions: 10 x 7.5 inches (standard) ─────────────────────────────
# officer default layout is 10 x 7.5 for widescreen

prs <- add_slide(prs, layout = "Blank", master = "Office Theme")

# ─── Colors ───────────────────────────────────────────────────────────────────
# We'll use ph_with for shapes. officer uses fp_* objects.

# Background – white slide is default

# Helper: RGB hex to officer color
col_navy    <- "#1B3A6B"
col_teal    <- "#0D7377"
col_slate   <- "#455A64"
col_seafoam <- "#00897B"
col_white   <- "#FFFFFF"
col_red     <- "#DC2626"
col_orange  <- "#D97706"
col_green   <- "#16A34A"
col_arrow   <- "#C0392B"
col_panel   <- "#EEF2F7"
col_legend  <- "#F7FAFC"
col_border  <- "#CBD5E0"

# ─── Layout constants (inches, origin top-left) ───────────────────────────────
LBL_X  <- 0.05
LBL_W  <- 1.25
BLX    <- 1.45    # box left x
BW     <- 2.85    # box width
CGAP   <- 0.5     # column gap
BRX    <- BLX + BW + CGAP   # 4.80
BTY    <- 1.15    # box top y
BH     <- 1.55    # box height
RGAP   <- 0.45    # row gap
BBY    <- BTY + BH + RGAP   # 3.15

# Slide dimensions
SW <- 10; SH <- 7.5

# ─── Add rectangle helper ─────────────────────────────────────────────────────
add_rect <- function(prs, x, y, w, h, fill, border = fill, border_w = 0) {
  ph_with(prs,
    value = block_list(
      fpar(ftext("", fp_text(font.size = 1)))
    ),
    location = ph_location(
      left = x, top = y, width = w, height = h,
      bg = fill,
      ln_color = border,
      ln_width = border_w
    )
  )
}

# officer's ph_location doesn't directly do shapes; we use add_shape via xml approach.
# Better: use officer's ph_with + external_img or gg_plot approach.
#
# Actually let me use officer's shape capabilities directly via add_shape.
# officer >= 0.4.0 has ph_location_type and we can add shapes via xml.
#
# The cleanest approach: use officer's slide_summary-based add approach.
# Let's use ph_with with location and block_list for text + separate shapes.
#
# Actually, the most straightforward way is to use the officer::ph_with
# with fpar/ftext for text boxes, and for colored rectangles we use
# officer's on_slide approach with xml_node for drawing ML shapes.
#
# Let me use a ggplot2-based approach instead: draw everything with ggplot2
# as a full-slide image, then embed it. This is simpler and gives full control.

# Close the pptx we started and use ggplot2
rm(prs)

library(ggplot2)
library(grid)
library(png)

# ─── Draw with ggplot2 ────────────────────────────────────────────────────────
# We'll create a plot (10 x 5.625 inch canvas = 16:9) and save as PNG,
# then embed into PPTX.

W <- 10; H <- 5.625

# Coordinate system: x in [0,10], y in [0, 5.625] with y increasing UPWARD
# (ggplot default). We'll flip mentally: top of slide = H, bottom = 0.

# Convert from "top-left inch coords" (ty = from top) to ggplot coords (y from bottom)
ty2y <- function(ty, h = 0) H - ty - h  # bottom y of a box at top-y=ty with height h

# Box positions in ggplot (xmin, xmax, ymin, ymax)
boxes <- list(
  list(label1="Smoked",   label2="Dark Filter",  fill=col_navy,
       xmin=BLX, xmax=BLX+BW, ymin=ty2y(BTY,BH), ymax=ty2y(BTY)),
  list(label1="Smoked",   label2="Light Filter", fill=col_teal,
       xmin=BRX, xmax=BRX+BW, ymin=ty2y(BTY,BH), ymax=ty2y(BTY)),
  list(label1="Unsmoked", label2="Dark Filter",  fill=col_slate,
       xmin=BLX, xmax=BLX+BW, ymin=ty2y(BBY,BH), ymax=ty2y(BBY)),
  list(label1="Unsmoked", label2="Light Filter", fill=col_seafoam,
       xmin=BRX, xmax=BRX+BW, ymin=ty2y(BBY,BH), ymax=ty2y(BBY))
)

p <- ggplot() +
  theme_void() +
  coord_cartesian(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# ── White background
p <- p + annotate("rect", xmin=0, xmax=W, ymin=0, ymax=H,
                  fill="white", color=NA)

# ── Title
p <- p + annotate("text",
  x = W/2, y = H - 0.35,
  label = "Cigarette Filter Compound Fate — Study Design",
  color = col_navy, size = 7.5, fontface = "bold", family = "sans", hjust = 0.5)

# ── Axis super-labels
mid_cols_x <- (BLX + BRX + BW) / 2
p <- p + annotate("text",
  x = mid_cols_x, y = H - 0.72,
  label = "IRRADIATION STATUS",
  color = "#999999", size = 3.3, fontface = "italic", hjust = 0.5)

# Column headers
p <- p + annotate("text",
  x = BLX + BW/2, y = H - 0.87,
  label = "DARK", color = col_navy, size = 4.5, fontface = "bold", hjust = 0.5)
p <- p + annotate("text",
  x = BRX + BW/2, y = H - 0.87,
  label = "LIGHT / IRRADIATED", color = col_teal, size = 4.5, fontface = "bold", hjust = 0.5)

# Row labels
mid_rows_y <- (ty2y(BTY) + ty2y(BBY, BH)) / 2 + BH/2
# Actually compute center of entire grid vertically
grid_top    <- ty2y(BTY)       # ggplot y of top of top boxes
grid_bottom <- ty2y(BBY, BH)   # ggplot y of bottom of bottom boxes
mid_grid_y  <- (grid_top + grid_bottom) / 2

p <- p + annotate("text",
  x = LBL_X + LBL_W/2, y = mid_grid_y,
  label = "SMOKING\nSTATUS",
  color = "#999999", size = 3.0, fontface = "italic", hjust = 0.5, lineheight = 0.9)

p <- p + annotate("text",
  x = LBL_X + LBL_W/2,
  y = (ty2y(BTY) + ty2y(BTY, BH)) / 2,
  label = "SMOKED", color = col_navy, size = 4.3, fontface = "bold", hjust = 0.5)

p <- p + annotate("text",
  x = LBL_X + LBL_W/2,
  y = (ty2y(BBY) + ty2y(BBY, BH)) / 2,
  label = "UNSMOKED", color = col_slate, size = 4.3, fontface = "bold", hjust = 0.5)

# ── Grid boxes
for (b in boxes) {
  p <- p + annotate("rect",
    xmin = b$xmin, xmax = b$xmax, ymin = b$ymin, ymax = b$ymax,
    fill = b$fill, color = NA)
  # Shadow effect: thin dark rect slightly offset (fake shadow)
  # Box label
  cy <- (b$ymin + b$ymax) / 2
  cx <- (b$xmin + b$xmax) / 2
  p <- p + annotate("text", x = cx, y = cy + 0.18,
    label = b$label1, color = "white", size = 5.8, fontface = "bold", hjust = 0.5)
  p <- p + annotate("text", x = cx, y = cy - 0.18,
    label = b$label2, color = "white", size = 4.2, hjust = 0.5)
}

# ── Comparison arrows
# Arrow gap centers in ggplot y coords
row_gap_center_y <- (ty2y(BTY, BH) + ty2y(BBY)) / 2  # center of vertical gap
# Actually: ty2y(BTY,BH) = bottom of top boxes, ty2y(BBY) = top of bottom boxes
top_box_bottom <- ty2y(BTY, BH)
bot_box_top    <- ty2y(BBY)
col_gap_right  <- BLX + BW   # right edge of left column
col_gap_left   <- BRX        # left edge of right column

# Arrow helper using geom_segment with arrow
arr_end <- arrow(length = unit(0.2, "cm"), type = "closed")
arr_begin <- arrow(length = unit(0.2, "cm"), type = "closed", ends = "first")

# ① Smoked Dark → Unsmoked Dark (vertical down, left col center)
c1x <- BLX + BW/2
p <- p + annotate("segment",
  x = c1x, xend = c1x,
  y = top_box_bottom - 0.02, yend = bot_box_top + 0.02,
  color = col_arrow, linewidth = 1.2,
  arrow = arrow(length = unit(0.18, "cm"), type = "closed", ends = "last"))
p <- p + annotate("text",
  x = c1x - 0.22, y = row_gap_center_y,
  label = "①", color = col_arrow, size = 3.8, fontface = "bold")

# ② Smoked Light → Unsmoked Light (vertical down, right col center)
c2x <- BRX + BW/2
p <- p + annotate("segment",
  x = c2x, xend = c2x,
  y = top_box_bottom - 0.02, yend = bot_box_top + 0.02,
  color = col_arrow, linewidth = 1.2,
  arrow = arrow(length = unit(0.18, "cm"), type = "closed", ends = "last"))
p <- p + annotate("text",
  x = c2x + 0.22, y = row_gap_center_y,
  label = "②", color = col_arrow, size = 3.8, fontface = "bold")

# ③ Smoked Light → Smoked Dark (horizontal left, top row)
c3y_top   <- ty2y(BTY)     # top of top boxes (ggplot)
c3y_bot   <- ty2y(BTY,BH)  # bottom of top boxes
c3y <- c3y_bot + (c3y_top - c3y_bot) * 0.65   # 65% up in top boxes
p <- p + annotate("segment",
  x = col_gap_right + 0.02, xend = col_gap_left - 0.02,
  y = c3y, yend = c3y,
  color = col_arrow, linewidth = 1.2,
  arrow = arrow(length = unit(0.18, "cm"), type = "closed", ends = "first"))
p <- p + annotate("text",
  x = (col_gap_right + col_gap_left)/2, y = c3y + 0.2,
  label = "③", color = col_arrow, size = 3.8, fontface = "bold")

# ④ Unsmoked Light → Unsmoked Dark (horizontal left, bottom row)
c4y_top   <- ty2y(BBY)
c4y_bot   <- ty2y(BBY,BH)
c4y <- c4y_bot + (c4y_top - c4y_bot) * 0.35   # 35% up in bottom boxes
p <- p + annotate("segment",
  x = col_gap_right + 0.02, xend = col_gap_left - 0.02,
  y = c4y, yend = c4y,
  color = col_arrow, linewidth = 1.2,
  arrow = arrow(length = unit(0.18, "cm"), type = "closed", ends = "first"))
p <- p + annotate("text",
  x = (col_gap_right + col_gap_left)/2, y = c4y - 0.2,
  label = "④", color = col_arrow, size = 3.8, fontface = "bold")

# ── Side panel: Comparison key
PX <- BRX + BW + 0.18
PW <- W - PX - 0.12
PY_top <- ty2y(BTY)       # top of panel (ggplot)
PY_bot <- ty2y(BBY, BH)   # bottom of panel

p <- p + annotate("rect",
  xmin = PX, xmax = PX + PW,
  ymin = PY_bot, ymax = PY_top + 0.05,
  fill = col_panel, color = col_border, linewidth = 0.4)

p <- p + annotate("text",
  x = PX + PW/2, y = PY_top - 0.05,
  label = "Comparisons", color = col_navy, size = 3.5, fontface = "bold", hjust = 0.5)

comp_data <- list(
  list(n="①", a="Smoked Dark",    b="Unsmoked Dark"),
  list(n="②", a="Smoked Light",   b="Unsmoked Light"),
  list(n="③", a="Smoked Light",   b="Smoked Dark"),
  list(n="④", a="Unsmoked Light", b="Unsmoked Dark")
)

panel_h <- PY_top - PY_bot - 0.3
item_h  <- panel_h / 4

for (i in seq_along(comp_data)) {
  item_y <- PY_top - 0.3 - (i - 0.5) * item_h
  cd <- comp_data[[i]]
  p <- p + annotate("text",
    x = PX + 0.2, y = item_y,
    label = cd$n, color = col_arrow, size = 4.2, fontface = "bold", hjust = 0.5)
  p <- p + annotate("text",
    x = PX + 0.38, y = item_y + item_h * 0.18,
    label = cd$a, color = "#1B3A6B", size = 2.7, fontface = "bold", hjust = 0)
  p <- p + annotate("text",
    x = PX + 0.38, y = item_y - item_h * 0.18,
    label = paste0("→ ", cd$b), color = "#374151", size = 2.6, hjust = 0)
}

# ── Legend bar at bottom
LY_top <- 0.95
LY_bot <- 0.12
p <- p + annotate("rect",
  xmin = 0.1, xmax = W - 0.1,
  ymin = LY_bot, ymax = LY_top,
  fill = col_legend, color = "#E2E8F0", linewidth = 0.4)

p <- p + annotate("text",
  x = 0.22, y = (LY_top + LY_bot)/2,
  label = "Per-Compound\nOutcome:", color = col_navy,
  size = 3.2, fontface = "bold", hjust = 0, lineheight = 0.9)

legend_items <- list(
  list(label="Persisted",   color=col_green,  desc="Compound remains present"),
  list(label="Transformed", color=col_orange, desc="Chemically changed form"),
  list(label="Removed",     color=col_red,    desc="No longer detected")
)
for (i in seq_along(legend_items)) {
  li <- legend_items[[i]]
  lx <- 1.9 + (i - 1) * 2.65
  ly <- (LY_top + LY_bot) / 2
  # Colored dot
  p <- p + annotate("point",
    x = lx + 0.12, y = ly,
    color = li$color, size = 4, shape = 16)
  # Label text
  p <- p + annotate("text",
    x = lx + 0.3, y = ly,
    label = paste0(li$label, "  —  ", li$desc),
    color = "#1F2937", size = 3.1, hjust = 0, fontface = "plain")
}

# ── Save as PNG then embed in PPTX ─────────────────────────────────────────────
tmp_png <- tempfile(fileext = ".png")
ggsave(tmp_png, plot = p, width = W, height = H, dpi = 200, bg = "white")

# Create PPTX with the plot image filling the slide
prs2 <- read_pptx()
prs2 <- add_slide(prs2, layout = "Blank", master = "Office Theme")
prs2 <- ph_with(prs2,
  value = external_img(tmp_png, width = W, height = H),
  location = ph_location(left = 0, top = 0, width = W, height = H))

print(prs2, target = out_path)
cat("Saved:", out_path, "\n")
