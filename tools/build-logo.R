## Regenerate the diffcp hex sticker artwork from the canonical SVG.
##
## The source-of-truth is man/figures/logo.svg.  This script renders
## the standard PNG variants (`logo.png` at 240 px wide for navbars,
## `logo@2x.png` at 480 px for Retina displays) and refreshes the
## pkgdown favicon set so the rendered site picks up any logo edits.
##
## Run from the package root:
##   Rscript tools/build-logo.R
##
## The `tools/` directory is .Rbuildignore'd, so this script is not
## shipped in the package tarball.

stopifnot(file.exists("man/figures/logo.svg"))
if (!requireNamespace("rsvg", quietly = TRUE)) {
  stop("Install the 'rsvg' R package: install.packages('rsvg').")
}

cat("Rendering man/figures/logo.png (240 px)...\n")
rsvg::rsvg_png("man/figures/logo.svg", "man/figures/logo.png",  width = 240)

cat("Rendering man/figures/logo@2x.png (480 px Retina)...\n")
rsvg::rsvg_png("man/figures/logo.svg", "man/figures/logo@2x.png", width = 480)

if (requireNamespace("pkgdown", quietly = TRUE)) {
  cat("Rebuilding pkgdown favicons from the new logo...\n")
  pkgdown::build_favicons(overwrite = TRUE)
} else {
  cat("Skipping pkgdown favicons (pkgdown not installed).\n")
}

## Color palette + design notes (kept here so anyone editing the
## SVG knows the canonical values):
##
##   Hex border .................. #5a5d62  (medium-dark slate)
##   Hex fill .................... #ffffff  (with soft drop shadow)
##   Cone outline / "data" dot ... #001d55  (navy; class .navy)
##   Cone fill (gradient top) .... #ffffff @ 0.96
##   Cone fill (gradient bottom).. #dcecff @ 0.86
##   Inner cone edges ............ #2f66aa  (medium blue; class .blue)
##   Inner cone faint edges ...... #9bb9df  (pale blue; class .pblue)
##   Plane fill .................. #eaf3ff @ 0.55
##   "perturbed" / accent ........ #6733a4  (violet; class .purple)
##
##   Wordmark ............... Avenir > Montserrat > Helvetica > Arial,
##                            font-weight 700, font-size 155, letter-
##                            spacing 9, navy->violet linearGradient
##                            split at the f|c boundary.
##   Subtitle ............... same sans family, 500/27/9-tracked, all-
##                            caps, accent rule between wordmark and
##                            subtitle (#6733a4 stroke-width 5).
##   Annotation labels ...... Georgia (serif) italic 34 px.
##
## Geometry:
##   - Hex polygon points "600,50 1076,325 1076,875 600,1150 124,875 124,325"
##     (pointy-top, 1200x1200 viewBox).
##   - Convex set is a second-order (Lorentz) cone
##     { (x,y,t) : sqrt(x^2 + y^2) <= t }: apex at (600, 720),
##     elliptical rim centered at (600, 240) with rx=260, ry=64.
##     The front rim arc is solid navy; the back rim arc is dashed
##     pale-blue to suggest the hidden far side of the cone.
##   - Affine plane is a tilted parallelogram (273,582)-(421,460)-
##     (899,460)-(779,582) slicing through the cone.  Cone ∩ plane
##     is drawn as the small ellipse cx=600, cy=520, rx=118, ry=22
##     -- the intersection of two convex sets is convex.
##   - Solid navy "data" dot at the optimum (600, 542) on the
##     cone-plane intersection; violet "perturbed" dot up the cone
##     at (700, 380); dashed violet bezier between them is the
##     (dx, dy, ds) direction returned by `derivative()`.
##   - "perturbed" label sits in the upper-right quadrant; the
##     pointer uses dasharray "32 22" for two long dashes + arrowhead.
