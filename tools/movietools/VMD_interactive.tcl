# =======================
# VMD_interactive.tcl
# =======================

# ---- Config ----
set START        1
set END          160

set ISO1         0.002
set ISO2         1.0

# Subsample volumetric grid used by Isosurface (0=full, 1..2 faster)
set STEP         0

# Isosurface transparency (0.0 transparent .. 1.0 opaque)
set ISO_OPACITY  0.5
set ISO_MATERIAL "Transparent"   ;# will be given ISO_OPACITY

# Colors for the four isosurfaces, in order: +ISO1, -ISO1, +ISO2, -ISO2
# Use names (e.g., red, blue, orange, cyan) or ColorIDs (ints).
set ISO_COLORS {red blue orange cyan}

# Show atoms as CPK background?
set DRAW_ATOMS   1

# Autoplay settings
set AUTOPLAY     1
set PLAY_STYLE   "Loop"     ;# Loop | Rock | Once
set PLAY_SPEED   0.25

# Optional: spin camera while playing
set SPIN_CAMERA  0
set SPIN_AXIS    y
set SPIN_STEP    0.5
# ----------------


# ---- Helpers ----
proc resolve_color {c} {
  if {[string is integer -strict $c]} {
    return $c
  } else {
    return [colorinfo index $c]
  }
}

# Lists of rep indices and their fixed isovalues (same order)
set isoRepIdxs {}
set isoRepVals {}

proc add_iso {val color_spec} {
  # Create rep
  mol representation Isosurface $val 0 $::STEP 0
  mol Material $::ISO_MATERIAL
  mol addrep top
  set rid [expr {[molinfo top get numreps] - 1}]

  # Lock color per-rep
  set cid [resolve_color $color_spec]
  mol modcolor $rid top ColorID $cid

  # Track this rep and its isovalue
  lappend ::isoRepIdxs $rid
  lappend ::isoRepVals $val
}

# ---- Load cubes as frames of one molecule ----
set updmol [mol new "out${START}.cube" type cube waitfor all]
for {set i [expr {$START + 1}]} {$i <= $END} {incr i} {
  mol addfile "out${i}.cube" type cube waitfor all
}

display rendermode GLSL     ;# alpha-blended transparency (if supported)
display depthsort on        ;# correct draw order for transparent reps
display antialias on        ;# optional: nicer edges


# ---- Optional atoms rep ----
mol delrep 0 top
if {$DRAW_ATOMS} {
  mol representation CPK 0.4 0.2 20 12
  mol color Name
  mol addrep top
}

# ---- Set material opacity used by iso reps ----
material change opacity $ISO_MATERIAL $ISO_OPACITY

# ---- Build four iso reps with explicit colors ----
set isoVals [list $ISO1 [expr {-$ISO1}] $ISO2 [expr {-$ISO2}]]
for {set k 0} {$k < 4} {incr k} {
  add_iso [lindex $isoVals $k] [lindex $ISO_COLORS $k]
}

# ---- Camera like your template ----
rotate x by -90
rotate y by -90
scale by 0.45

# ---- Keep iso reps synced to current frame WITHOUT touching color ----
proc update_iso {args} {
  set f [molinfo $::updmol get frame]
  for {set j 0} {$j < [llength $::isoRepIdxs]} {incr j} {
    set rid  [lindex $::isoRepIdxs $j]
    set isov [lindex $::isoRepVals $j]
    # Only change the style params (isovalue, volumeID/frame, step); color stays put
    mol modstyle $rid $::updmol Isosurface $isov $f $::STEP 0
  }
}
trace variable vmd_frame($updmol) w update_iso

# ---- Autoplay ----
if {$AUTOPLAY} {
  animate style $PLAY_STYLE
  animate speed $PLAY_SPEED
  animate goto $START
  display update on
  animate forward
}

# ---- Optional camera spin ----
proc cam_spin {{axis y} {step 0.5}} {
  rotate $axis by $step
  after 33 [list cam_spin $axis $step]
}
if {$SPIN_CAMERA} {
  cam_spin $SPIN_AXIS $SPIN_STEP
}

# --- Notes (Tk Console):
#   mol modcolor 1 top ColorID [colorinfo index purple]
#   material change opacity $ISO_MATERIAL 0.3
