# source me in VMD terminal

set updmol [mol new {outTSTARTKEY.cube} type cube waitfor all]
for { set i TSTARTP1KEY } { $i <= TENDKEY } { incr i } {
	 set myfile out$i.cube 
	 mol addfile $myfile type cube waitfor all 
}

mol delrep 0 top
mol representation  CPK 0.4 0.2 20 12
mol color Name
mol addrep top

mol representation Isosurface  ISOKEY1 0.0 0.0 0.0
mol Material Transparent
mol color ColorID 21
mol addrep top

mol representation Isosurface -ISOKEY1 0.0 0.0 0.0
mol Material Transparent
mol color ColorID 1
mol addrep top


mol representation Isosurface  ISOKEY2 0.0 0.0 0.0
mol Material Transparent
mol color ColorID 3
mol addrep top

mol representation Isosurface -ISOKEY2 0.0 0.0 0.0
mol Material Transparent
mol color ColorID 1
mol addrep top





set updrep1 [mol repname top 1]
set updrep2 [mol repname top 2]
set updrep3 [mol repname top 3]
set updrep4 [mol repname top 4]


# Camera Adjustments
# Rotate to look down y-axis
rotate x by -90
rotate y by -90
# Zoom out
scale by 0.45

proc update_iso {args} {
   global updmol
   global updrep1
   global updrep2
   global updrep3
   global updrep4

   set repid1 [mol repindex $updmol $updrep1]
   if { $repid1 < 0} { return }

   set frame [molinfo $updmol get frame]

   lassign [molinfo $updmol get "{rep $repid1}"] rep
   mol representation [lreplace $rep 2 2 $frame]
   mol color ColorID 21
   mol modrep $repid1 $updmol

   set repid2 [mol repindex $updmol $updrep2]
   if { $repid2 < 0} { return }

   set frame [molinfo $updmol get frame]

   lassign [molinfo $updmol get "{rep $repid2}"] rep
   mol representation [lreplace $rep 2 2 $frame]
   mol color ColorID 1
   mol modrep $repid2 $updmol


   set repid3 [mol repindex $updmol $updrep3]
   if { $repid3 < 0} { return }

   set frame [molinfo $updmol get frame]

   lassign [molinfo $updmol get "{rep $repid3}"] rep
   mol representation [lreplace $rep 2 2 $frame]
   mol color ColorID 21
   mol modrep $repid3 $updmol

   set repid4 [mol repindex $updmol $updrep4]
   if { $repid4 < 0} { return }

   set frame [molinfo $updmol get frame]

   lassign [molinfo $updmol get "{rep $repid4}"] rep
   mol representation [lreplace $rep 2 2 $frame]
   mol color ColorID 1
   mol modrep $repid4 $updmol

}

trace variable vmd_frame($updmol) w update_iso

set renderwidth 800
set renderheight 600

for {set i TSTARTKEY} {$i <= TENDKEY} {incr i} {
   animate goto $i
   set filename "render$i.bmp"
   render TachyonInternal $filename
}



quit
